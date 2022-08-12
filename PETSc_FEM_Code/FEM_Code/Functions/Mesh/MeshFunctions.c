#ifndef MeshFunctions_H
#define MeshFunctions_H


PetscErrorCode MeshFindGhostNodes(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements,
  const PetscInt N_U_total_glo, PetscInt *OWNER_RANK, PetscInt *OWNER_LI, PetscInt *OWNER_T)
{
  PetscErrorCode ierr;

  PetscInt N_recv;
  PetscInt *recv_GI, *recv_t;

  PetscMPIInt size = system->size;
  PetscMPIInt rank = system->rank;     
  NodesCountU(nodes, rank);
  if (size >1)
  {
    for(PetscInt i=0; i<nodes->N_U_total_loc; i++) 
      if(nodes->t[i] >1){
        SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in FindGhostNodes"); }

    for(PetscInt i=0; i<N_U_total_glo; i++) 
      OWNER_RANK[i] = -1, OWNER_LI[i] = -1, OWNER_T[i]=0; 

    for(PetscInt r=0; r<size; r++)
    {
      if(rank == r)
      {
        ierr = MPI_Bcast(&nodes->N_U_total_loc, 1,          MPIU_INT, r, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
        ierr = MPI_Bcast(nodes->GI_U, nodes->N_U_total_loc, MPIU_INT, r, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
        ierr = MPI_Bcast(nodes->t,    nodes->N_U_total_loc, MPIU_INT, r, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 

        for(PetscInt i=0; i<nodes->N_U_total_loc; i++) {
          const PetscInt gi = nodes->GI_U[i];
          if (OWNER_RANK[gi] == -1)
            OWNER_RANK[gi] = r, OWNER_LI[gi] = i; 
          OWNER_T[gi] += nodes->t[i];
        }
      }
      else
      {
        ierr = MPI_Bcast(&N_recv, 1,        MPIU_INT, r, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
        ierr = PetscMalloc2(N_recv, &recv_GI, N_recv, &recv_t); CHKERRQ(ierr);
        ierr = MPI_Bcast(recv_GI, N_recv, MPIU_INT, r, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
        ierr = MPI_Bcast(recv_t,  N_recv, MPIU_INT, r, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 

        for(PetscInt i=0; i<N_recv; i++) {
          const PetscInt gi = recv_GI[i];
          if (OWNER_RANK[gi] == -1)
            OWNER_RANK[gi] = r, OWNER_LI[gi] = i;
          OWNER_T[gi] += recv_t[i];
        }

        ierr = PetscFree2(recv_GI, recv_t); CHKERRQ(ierr);
      }
    }

    for(PetscInt i=0; i<N_U_total_glo; i++)
    {
      if (OWNER_RANK[i] == -1 || OWNER_LI[i] == -1) {
        SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in FindGhostNodes"); }
        OWNER_T[i] = PetscMin(1,OWNER_T[i]); 
    }

    for(PetscInt i=0; i<nodes->N_U_total_loc; i++){
      const PetscInt gi = nodes->GI_U[i];
      if (OWNER_RANK[gi] != rank)
        nodes->t[i] = OWNER_T[gi] + 2;
      else
        nodes->t[i] = OWNER_T[gi];
    }
    NodesCountU(nodes, rank);
  }

  return 0;
}

PetscErrorCode MeshSaveUBoundaryAndGhostLI(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements, 
  const PetscInt N_U_total_glo, PetscInt *OWNER_RANK, PetscInt *OWNER_LI)
{
  PetscInt i, j, k;

  NodesCountU(nodes, system->rank);
  NodesMallocUBoundaryArrays(nodes, system->rank);
  NodesMallocUGhostArrays(nodes, system->rank);

  j = 0, k = 0;
  for(i = 0; i < nodes->N_U_total_loc; ++i)
    if(nodes->t[i] == 1)                           
      nodes->LI_U_bound[j] = i, j++;
    else if (nodes->t[i] == 2 || nodes->t[i] == 3) 
      nodes->LI_U_ghost[k] = i, k++;

  if (j != nodes->N_U_bound_loc){
    PetscErrorPrintf("The number of boundary nodes found is %d it should be %d\n", j, nodes->N_U_bound_loc); 
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in SaveUBoundaryAndGhostLI"); }
  if (k != nodes->N_U_ghost_loc){
    PetscErrorPrintf("The number of ghost nodes found is %d it should be %d\n", k, nodes->N_U_ghost_loc); 
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in SaveUBoundaryAndGhostLI"); }
  
  for(i=0; i < elements->N_bound_loc; i++)
    for (j=0; j < elements->U_S_bound; j++)
    {
      for(k=0; k < nodes->N_U_bound_loc; k++)
        if(nodes->LI_U_bound[k] == elements->bound[i][j])
          break;
      if (k != nodes->N_U_bound_loc) 
        nodes->Gr_U_bound[k] = elements->Gr_bound[i];
    }
  
  for(i = 0; i < nodes->N_U_ghost_loc; ++i)
  {
    const PetscInt li = nodes->LI_U_ghost[i];
    const PetscInt gi = nodes->GI_U[li];
    nodes->Ra_U_ghost[i]    = OWNER_RANK[gi];
    nodes->Ra_LI_U_ghost[i] = OWNER_LI[gi];
  }

  PetscPrintf(PETSC_COMM_WORLD,"  Done saving the U boundary and ghost LI arrays on all ranks\n");
  return 0;
}



PetscErrorCode MeshCreatePNodesFromUNodes(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements)
{
  PetscInt e, i, j, k, LI_P;

  nodes->N_P_total_loc = 0;
  nodes->N_P_bound_loc = 0;
  nodes->N_P_inter_loc = 0;
  nodes->N_P_ghost_loc = 0;
  
  for(i = 0; i < nodes->N_U_total_loc; ++i)
    nodes->LI_UtoP[i] = -1;
  
  for(e=0; e < elements->N_inter_loc; e++)
    for(i=0; i < elements->P_S_inter; i++)
      nodes->LI_UtoP[elements->inter[e][i]] = 0; 
  
  for(e=0; e < elements->N_bound_loc; e++)
    for(i=0; i < elements->P_S_bound; i++)
      nodes->LI_UtoP[elements->bound[e][i]] = 1; 

  LI_P = 0; 
  for(i = 0; i < nodes->N_U_total_loc; ++i)
    if (nodes->LI_UtoP[i] != -1)
    {
      nodes->N_P_total_loc++;
      if (nodes->t[i] >= 2)
        nodes->N_P_ghost_loc++;
      else if(nodes->t[i] == 1)
        nodes->N_P_bound_loc++;
      else if (nodes->t[i] == 0)
        nodes->N_P_inter_loc++;
      nodes->LI_UtoP[i] = LI_P;
      LI_P++;
    }

  NodesMallocPArraysForAll(nodes, system->rank);
  NodesMallocPGhostArrays(nodes, system->rank);

  k = 0; 
  for(i = 0; i < nodes->N_U_total_loc; ++i)
    if (nodes->LI_UtoP[i] != -1)
    {
      LI_P = nodes->LI_UtoP[i];
      if (nodes->t[i] >= 2)
      {
        nodes->LI_P_ghost[k] = LI_P; 
        for(j =0; j< nodes->N_U_ghost_loc; j++)
          if(nodes->LI_U_ghost[j] == i)
            break;
        if(j==nodes->N_U_ghost_loc){
          SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in CreatePNodesFromUNodes"); }
        nodes->Ra_P_ghost[k]      = nodes->Ra_U_ghost[j];
        nodes->Ra_ULI_P_ghost[k]  = nodes->Ra_LI_U_ghost[j];
        k++; 
      }
      nodes->LI_PtoU[LI_P] = i;          
      nodes->GI_P[LI_P] = nodes->GI_U[i];
    }

  PetscPrintf(PETSC_COMM_WORLD,"  Done creating the P nodes on all ranks\n");
  return 0;
}


PetscErrorCode MeshReOrganizeNodesGI(System_Struct *system, Nodes_Struct *nodes)
{
  PetscErrorCode ierr;
  PetscInt i, r, gi_U, gi_P;
  PetscInt *N_total_glo;  
  PetscMPIInt size = system->size;
  PetscMPIInt rank = system->rank;        
  nodes->N_cores[rank][0] = nodes->N_P_total_loc;
  nodes->N_cores[rank][1] = nodes->N_P_ghost_loc;
  nodes->N_cores[rank][2] = nodes->N_U_total_loc;  
  nodes->N_cores[rank][3] = nodes->N_U_ghost_loc;

  for (PetscInt r=0; r<size; r++){
    ierr = MPI_Bcast(nodes->N_cores[r], 4, MPIU_INT, r, PETSC_COMM_WORLD ); CHKERRMPI(ierr);}
  
  nodes->GI_U_start = 0, nodes->GI_P_start = 0;
  for (r=0; r < rank; r++)
    nodes->GI_U_start += nodes->N_cores[r][2] - nodes->N_cores[r][3],
    nodes->GI_P_start += nodes->N_cores[r][0] - nodes->N_cores[r][1];

  gi_U = nodes->GI_U_start;
  for(i=0; i<nodes->N_U_total_loc; i++)
  {
    if (nodes->t[i] < 2) 
    {
      nodes->GI_U[i] = gi_U;         
      gi_U++;            
    }
  }

  gi_P = nodes->GI_P_start;
  for(i=0; i<nodes->N_P_total_loc; i++)
  {
    if (nodes->t[nodes->LI_PtoU[i]] < 2) 
    {
      nodes->GI_P[i] = gi_P;         
      gi_P++;            
    }
  }

  if (size == 1)
  {
    nodes->N_P_total_glo = nodes->N_P_total_loc;
    nodes->N_P_inter_glo = nodes->N_P_inter_loc;
    nodes->N_P_bound_glo = nodes->N_P_bound_loc;
    nodes->N_U_total_glo = nodes->N_U_total_loc;
    nodes->N_U_inter_glo = nodes->N_U_inter_loc;
    nodes->N_U_bound_glo = nodes->N_U_bound_loc;
  }
  else
  {
    ierr = PetscMalloc1(2, &N_total_glo); CHKERRQ(ierr);
    if (rank == size-1) 
    {
      for(i=nodes->N_P_total_loc-1; i>=0; i--)
        if(nodes->t[nodes->LI_PtoU[i]]<2)
          break;
      nodes->N_P_total_glo = nodes->GI_P[i] + 1;

      for(i=nodes->N_U_total_loc-1; i>=0; i--)
        if(nodes->t[i]<2)
          break;
      nodes->N_U_total_glo = nodes->GI_U[i] + 1;

      N_total_glo[0] = nodes->N_P_total_glo, N_total_glo[1] = nodes->N_U_total_glo;
      ierr = MPI_Bcast(N_total_glo, 2, MPIU_INT, size-1, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    }
    else
    {
      ierr = MPI_Bcast(N_total_glo, 2, MPIU_INT, size-1, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      nodes->N_P_total_glo = N_total_glo[0], nodes->N_U_total_glo = N_total_glo[1];
    }
    ierr = PetscFree(N_total_glo); CHKERRQ(ierr);
    SumIntAcrossRanks(size, rank, &nodes->N_P_inter_glo, nodes->N_P_inter_loc);
    SumIntAcrossRanks(size, rank, &nodes->N_P_bound_glo, nodes->N_P_bound_loc);
    SumIntAcrossRanks(size, rank, &nodes->N_U_inter_glo, nodes->N_U_inter_loc);
    SumIntAcrossRanks(size, rank, &nodes->N_U_bound_glo, nodes->N_U_bound_loc);

    if (nodes->N_P_inter_glo + nodes->N_P_bound_glo != nodes->N_P_total_glo){
      PetscErrorPrintf("Error: nodes->N_P_inter_glo + nodes->N_P_bound_glo (%d) != nodes->N_P_total_glo (%d) on rank %d\n",
       nodes->N_P_inter_glo + nodes->N_P_bound_glo, nodes->N_P_total_glo, rank); 
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in ReOrganizeNodesGI"); }

    if (nodes->N_U_inter_glo + nodes->N_U_bound_glo != nodes->N_U_total_glo){
      PetscErrorPrintf("Error: nodes->N_U_inter_glo + nodes->N_U_bound_glo (%d) != nodes->N_U_total_glo (%d) on rank %d\n",
       nodes->N_U_inter_glo + nodes->N_U_bound_glo, nodes->N_U_total_glo, rank); 
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in ReOrganizeNodesGI"); }
  }

  PetscPrintf(PETSC_COMM_WORLD,"  Done reorganizing the node global indices on all ranks\n");  
  return 0;
}





PetscErrorCode MeshShareNewGhostGI(System_Struct *system, Nodes_Struct *nodes)
{
  PetscErrorCode ierr;
  PetscInt s, r, i;
  PetscInt *Ra_U_send, *Ra_P_send, *I1, *I2;
  PetscInt **Ra_LI_U_send, **Ra_ULI_P_send;  
  PetscInt **LI_U_ghost, **LI_P_ghost;      

  PetscInt *LI_U_receive, *ULI_P_receive;
  PetscInt *GI_U_send, *GI_P_send, *GI_P_receive, *GI_U_receive;
  PetscInt Ra_U_receive, Ra_P_receive;

  PetscMPIInt size = system->size;
  PetscMPIInt rank = system->rank;       
  MPI_Request Request;
  MPI_Status Status;

  ierr = PetscMalloc4(rank, &Ra_U_send, rank ,&Ra_P_send, rank, &I1, rank ,&I2); CHKERRQ(ierr);
  for(s = 0; s < rank; ++s)
    Ra_U_send[s] = 0, Ra_P_send[s] = 0, I1[s] = 0, I2[s] = 0;

  for (i=0; i < nodes->N_P_ghost_loc; i++)  Ra_P_send[nodes->Ra_P_ghost[i]]++;
  for (i=0; i < nodes->N_U_ghost_loc; i++)  Ra_U_send[nodes->Ra_U_ghost[i]]++; 
  
  ierr = PetscMalloc4(rank, &Ra_LI_U_send, rank, &Ra_ULI_P_send, rank, &LI_U_ghost, rank, &LI_P_ghost);    CHKERRQ(ierr);
  for(s = 0; s < rank; ++s){
    ierr = PetscMalloc4( Ra_U_send[s], &Ra_LI_U_send[s], Ra_P_send[s], &Ra_ULI_P_send[s], 
      Ra_U_send[s], &LI_U_ghost[s], Ra_P_send[s], &LI_P_ghost[s]); CHKERRQ(ierr); }
  
  for (i=0; i < nodes->N_P_ghost_loc; i++) { 
    const PetscInt Ra = nodes->Ra_P_ghost[i];
    Ra_ULI_P_send[Ra][I1[Ra]] = nodes->Ra_ULI_P_ghost[i];
    LI_P_ghost[Ra][I1[Ra]] = nodes->LI_P_ghost[i];
    I1[Ra]++; 
    }

  for (i=0; i < nodes->N_U_ghost_loc; i++) { 
    const PetscInt Ra = nodes->Ra_U_ghost[i];
    Ra_LI_U_send[Ra][I2[Ra]] = nodes->Ra_LI_U_ghost[i];
    LI_U_ghost[Ra][I2[Ra]] = nodes->LI_U_ghost[i];
    I2[Ra]++; }

  for (s=1; s < size; s++)
  {
    if (rank == s) 
    {
      
      for (r=0; r < s; ++r)  
      {
        ierr = MPI_Isend( &Ra_P_send[r], 1, MPIU_INT, r, 1000+s, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
        ierr = MPI_Isend( Ra_ULI_P_send[r], Ra_P_send[r], MPIU_INT, r, 2000+s, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
        ierr = MPI_Isend( &Ra_U_send[r], 1, MPIU_INT, r, 3000+s, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
        ierr = MPI_Isend( Ra_LI_U_send[r], Ra_U_send[r], MPIU_INT, r, 4000+s, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request);   CHKERRMPI(ierr);
      }
      
      
      for (r=0; r < s; ++r)  
      {
        ierr = PetscMalloc2(Ra_P_send[r], &GI_P_receive, Ra_U_send[r], &GI_U_receive); CHKERRMPI(ierr);
        ierr = MPI_Irecv(GI_P_receive, Ra_P_send[r], MPIU_INT, r, 5000+s, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
        ierr = MPI_Irecv(GI_U_receive, Ra_U_send[r], MPIU_INT, r, 6000+s, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
        for(i=0; i< Ra_P_send[r]; i++) nodes->GI_P[LI_P_ghost[r][i]] = GI_P_receive[i];
        for(i=0; i< Ra_U_send[r]; i++) nodes->GI_U[LI_U_ghost[r][i]] = GI_U_receive[i];
        ierr = PetscFree2(GI_P_receive, GI_U_receive); CHKERRQ(ierr);
      }
    }
    else if (rank < s)
    {
      
      ierr = MPI_Irecv(&Ra_P_receive, 1, MPIU_INT, s, 1000+s, PETSC_COMM_WORLD, &Request);            CHKERRMPI(ierr); 
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request);           CHKERRMPI(ierr);
      ierr = PetscMalloc2(Ra_P_receive, &ULI_P_receive, Ra_P_receive, &GI_P_send); 
      ierr = MPI_Irecv(ULI_P_receive, Ra_P_receive, MPIU_INT, s, 2000+s, PETSC_COMM_WORLD, &Request);  CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request);           CHKERRMPI(ierr);
      ierr = MPI_Irecv(&Ra_U_receive, 1, MPIU_INT, s, 3000+s, PETSC_COMM_WORLD, &Request);            CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request);           CHKERRMPI(ierr);
      ierr = PetscMalloc2(Ra_U_receive, &LI_U_receive, Ra_U_receive, &GI_U_send); 
      ierr = MPI_Irecv(LI_U_receive, Ra_U_receive, MPIU_INT, s, 4000+s, PETSC_COMM_WORLD, &Request);  CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request);           CHKERRMPI(ierr);

      for(i=0; i<Ra_P_receive; i++) GI_P_send[i] = nodes->GI_P[nodes->LI_UtoP[ULI_P_receive[i]]];
      for(i=0; i<Ra_U_receive; i++) GI_U_send[i] = nodes->GI_U[LI_U_receive[i]];   
        
      ierr = MPI_Isend( GI_P_send, Ra_P_receive, MPIU_INT, s, 5000+s, PETSC_COMM_WORLD, &Request);  CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request);         CHKERRMPI(ierr);
      ierr = MPI_Isend( GI_U_send, Ra_U_receive, MPIU_INT, s, 6000+s, PETSC_COMM_WORLD, &Request);  CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request);         CHKERRMPI(ierr);
      ierr = PetscFree2(ULI_P_receive, GI_P_send); CHKERRQ(ierr);
      ierr = PetscFree2(LI_U_receive, GI_U_send); CHKERRQ(ierr);
    }
    ierr = PetscBarrier(NULL); CHKERRQ(ierr);
  }

  
  for(s = 0; s < rank; ++s)  ierr = PetscFree4( Ra_LI_U_send[s], Ra_ULI_P_send[s], LI_U_ghost[s], LI_P_ghost[s]); CHKERRQ(ierr);
  ierr = PetscFree4(Ra_LI_U_send, Ra_ULI_P_send, LI_U_ghost, LI_P_ghost);  CHKERRQ(ierr);
  ierr = PetscFree4(Ra_U_send, Ra_P_send, I1, I2);   CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"    Done reorganizing the ghost node global indices on all ranks\n");  
  return 0;
}

PetscErrorCode MeshCreateP1bElements(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements)
{
 
  PetscInt e, n, i, N;
  PetscScalar xc, yc ,zc;
 
  PetscScalar *x_n, *y_n ,*z_n;
  PetscInt    *t_n, *GI_U_n, **e_n; 
  const PetscInt elem_size = elements->U_S_inter;

  PetscInt Old_N_U_total_loc = nodes->N_U_total_loc;
  PetscInt New_N_U_total_loc = nodes->N_U_total_loc + elements->N_inter_loc;
  
  PetscMalloc5(New_N_U_total_loc, &x_n, New_N_U_total_loc, &y_n,    New_N_U_total_loc, &z_n, 
               New_N_U_total_loc, &t_n, New_N_U_total_loc, &GI_U_n);
  PetscMalloc1(New_N_U_total_loc, &e_n);
  for(i=0; i<New_N_U_total_loc; i++) PetscMalloc1(3, &e_n[i]);  

  for(i=0; i<Old_N_U_total_loc; i++) 
  {
    e_n[i][0]=nodes->e[i][0], e_n[i][1]=nodes->e[i][1], e_n[i][2]=nodes->e[i][2];
    x_n[i] = nodes->x[i];
    y_n[i] = nodes->y[i];
    z_n[i] = nodes->z[i];
    t_n[i] = nodes->t[i];
    GI_U_n[i] = nodes->GI_U[i];
  }

  N = Old_N_U_total_loc; 
  for (e=0; e<elements->N_inter_loc; e++) 
  {
    
    xc=0, yc=0, zc=0;
    for(i=0; i<elem_size-1; i++ )     {
      n = elements->inter[e][i];
      xc += x_n[n]; yc += y_n[n]; zc += z_n[n];
    }
    xc = xc/(elem_size-1); yc = yc/(elem_size-1); zc = zc/(elem_size-1);
    
    x_n[N] = xc; y_n[N] = yc; z_n[N] = zc; t_n[N] = 0, e_n[N][0]=1, e_n[N][1]=0, e_n[N][2]=0;
    GI_U_n[N] = N;

    elements->inter[e][elem_size-1] = N;
    N++;
  }
 
  for(i=0; i<Old_N_U_total_loc; i++) PetscFree(nodes->e[i]);
  PetscFree(nodes->e);
  PetscFree6(nodes->x, nodes->y, nodes->z, nodes->t, nodes->GI_U, nodes->LI_UtoP);
  for(i=0; i<system->size; i++)  PetscFree(nodes->N_cores[i]);
  PetscFree(nodes->N_cores);
  
  nodes->N_U_total_loc = New_N_U_total_loc;
  NodesMallocUArraysForAll(nodes, system->size, system->rank);
  for(i=0; i<New_N_U_total_loc; i++) 
  {
    nodes->e[i][0] = e_n[i][0];
    nodes->e[i][1] = e_n[i][1];
    nodes->e[i][2] = e_n[i][2];
    nodes->x[i]    = x_n[i];
    nodes->y[i]    = y_n[i];
    nodes->z[i]    = z_n[i];
    nodes->t[i]    = t_n[i];

    nodes->GI_U[i] = GI_U_n[i];
  }

  PetscFree5(x_n, y_n, z_n, t_n, GI_U_n);
  for(i=0; i<New_N_U_total_loc; i++) PetscFree(e_n[i]);
  PetscFree(e_n);
  return 0;
}



PetscErrorCode MeshCheckIfCorrectMesh(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements,
  const PetscInt N_U_total_glo, const PetscInt *OWNER_RANK, const PetscInt *OWNER_LI)
{  
  if (N_U_total_glo + elements->N_inter_glo != nodes->N_U_total_glo){
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in MeshCheckIfCorrectMesh"); }
  
  if (N_U_total_glo != nodes->N_P_total_glo){
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in MeshCheckIfCorrectMesh"); }
  
  for(PetscInt ei=0; ei<elements->N_inter_loc; ei++)
    CheckBaseCompatibility(system, elements->U_elem_type, ei, nodes->x, nodes->y, nodes->z, elements->inter);

  for(PetscInt ei=0; ei<elements->N_inter_loc; ei++)
    CheckBaseCompatibility(system, elements->P_elem_type, ei, nodes->x, nodes->y, nodes->z, elements->inter);
  
  for(PetscInt eb=0; eb<elements->N_bound_loc; eb++)
    CheckBoundaryCompatibility(system, elements->U_elem_type, elements->N_bound_loc, eb, nodes->x, nodes->y, nodes->z, elements->bound);
  
  for(PetscInt eb=0; eb<elements->N_bound_loc; eb++)
    CheckBoundaryCompatibility(system, elements->U_elem_type, elements->N_bound_loc, eb, nodes->x, nodes->y, nodes->z, elements->bound);
  
  for(PetscInt i=0; i < nodes->N_U_ghost_loc; i++)
  {
    PetscInt r, GI_Largest;
    const PetscInt Li   = nodes->LI_U_ghost[i];
    const PetscInt Ra_U = nodes->Ra_U_ghost[i];
    const PetscInt GI_U = nodes->GI_U[Li];

    GI_Largest = 0;
    
    for(r=0; r<system->size; r++) {
      GI_Largest += nodes->N_cores[r][2] - nodes->N_cores[r][3];
      if (GI_U < GI_Largest)
        break; 
    }
    if (r>=system->size){
      PetscErrorPrintf("Check the ghost node Li_U=%d, has a GI_U=%d greater than %d\n", Li, GI_U, GI_Largest); 
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in MeshCheckIfCorrectMesh"); }
    else if (r != Ra_U){
      PetscErrorPrintf("Check the ghost node Li_U=%d, has an owner rank of Ra_U=%d which is different from the rank found %d\n", Li, Ra_U, r); 
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in MeshCheckIfCorrectMesh"); }
  } 

  for(PetscInt i=0; i < nodes->N_P_ghost_loc; i++)
  {
    PetscInt r, GI_Largest;
    const PetscInt Li   = nodes->LI_P_ghost[i];
    const PetscInt Ra_P = nodes->Ra_P_ghost[i];
    const PetscInt GI_P = nodes->GI_P[Li];

    GI_Largest = 0;
    
    for(r=0; r<system->size; r++) {
      GI_Largest += nodes->N_cores[r][0] - nodes->N_cores[r][1];
      if (GI_P < GI_Largest)
        break; 
    }
    if (r>=system->size){
      PetscErrorPrintf("Check the ghost node Li_P=%d, has a GI_P=%d greater than %d\n", Li, GI_P, GI_Largest); 
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in MeshCheckIfCorrectMesh"); }
    else if (r != Ra_P){
      PetscErrorPrintf("Check the ghost node Li_P=%d, has an owner rank of Ra_P=%d which is different from the rank found %d\n", Li, Ra_P, r); 
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in MeshCheckIfCorrectMesh"); }
  }

  PetscPrintf(PETSC_COMM_WORLD,"  Mesh has passed all the preliminary tests!\n");
  NodesPrintInfoMini(system, nodes, PETSC_TRUE);
  ElementsPrintInfo(system, elements, PETSC_TRUE);

  return 0;
}



#include "VTKFunctions.c"
#include "GmshFunctions.c"


PetscErrorCode CreateFluidDomain(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements)
{  

  if (system->FLMesh_Type == GMSH)
    GMSHCreateFluidMesh(system, nodes, elements);
  else if (system->FLMesh_Type == FOAM){ 
    PetscErrorPrintf("FOAM Meshes are not defined yet. The only fluid mesh formats that are defined are: GMSH (msh).\n");
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Mesh type unknown.");  }
  else { 
    PetscErrorPrintf("Unknown mesh format in input file. The only fluid mesh formats that are defined are: GMSH (msh).\n");
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Mesh type unknown.");  }

  NodesPrintInfo(system, nodes, print_nodes);
  ElementsPrintInfo(system, elements, print_elements);
  
  if (write_fluid_mesh) 
    VTKWriteFluidMesh(system, nodes, elements);

  return 0;
}

PetscErrorCode CreateSolidDomain(System_Struct *system, Solid_Struct *solid)
{  

  if (system->SDMesh_Type == VTK)
    VTKCreateSolidMesh(system, solid);
  else { 
    PetscErrorPrintf("Unknown mesh format in input file. The only solid mesh formats that are defined are: VTK (vtk).\n");
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Mesh type unknown.");  }

  SolidConvertTrs2Lns(system, solid);

  SolidPrintInfo(system, solid, print_solid);

  return 0;
}

#endif