#ifndef CreateFullGlobalVectorsAndMatrices_C
#define CreateFullGlobalVectorsAndMatrices_C

PetscErrorCode CreateGlobalFullVectorWithGhostCells( Vec *V, const PetscInt glo, const PetscInt loc, 
  const PetscInt ghost_loc, const PetscInt *GI_V_ghost, const PetscMPIInt size)
{
  PetscErrorCode ierr;
  if(size > 1)
    {ierr = VecCreateGhost(PETSC_COMM_WORLD, loc, glo, ghost_loc, GI_V_ghost, V); CHKERRQ(ierr);}
  else
    {ierr = VecCreateSeq(PETSC_COMM_WORLD, loc, V); CHKERRQ(ierr);}

  ierr = VecSet(*V,0.0); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"  Finished creating the global vectors with global size %i on all ranks.\n", glo);
  
  return 0;
}

PetscErrorCode CalculateTheNumberOfConnections(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements,
  PetscInt **ln_M, PetscInt **ln_S, const PetscInt A_Width)
{
  PetscErrorCode ierr;

  const PetscInt size = system->size;
  const PetscInt rank = system->rank;

  PetscBool  Found = PETSC_FALSE;
  PetscInt i, j, k, u, r;

  PetscInt **gn_M, **gn_S; 
  PetscInt size_gn_M, size_gn_S; 

  PetscInt **rgn_M, **rgn_S; 
  PetscInt rsize_gn_M, rsize_gn_S; 

  PetscInt GI_start, GI_end; 

  PetscInt li;
  
  MPI_Request Request;
  MPI_Status Status;
  int flg;

  GI_start = 0;
  for (r=0; r < rank; r++)  
    GI_start += nodes->N_cores[r][2] - nodes->N_cores[r][3];
  GI_end = GI_start + nodes->N_cores[rank][2] - nodes->N_cores[rank][3];

  PetscPrintf(PETSC_COMM_WORLD,"  Countting the number of connections on all ranks.\n");
 
  for (i=0; i<elements->N_inter_loc; i++) 
  {
    for (j=0; j<elements->U_S_inter; j++) 
    {
      const PetscInt n = elements->inter[i][j];
      
      ln_M[n][0] = nodes->GI_U[n];
      ln_S[n][0] = nodes->GI_U[n];
    
      for(k=0; k<elements->P_S_inter; k++) 
      {
        const PetscInt m = elements->inter[i][k];
        const PetscInt gi_m = nodes->GI_U[m];
        
        Found = PETSC_FALSE;
        for(u=0; u < ln_M[n][1]; u++) 
          if(ln_M[n][u+2]==gi_m) {
            Found = PETSC_TRUE;
            break;}
        if (!Found) 
        { ln_M[n][1]++; ln_M[n][ln_M[n][1]+1] = gi_m; }
      }

      for(k=elements->P_S_inter; k<elements->U_S_inter; k++) 
      {
        const PetscInt m = elements->inter[i][k];
        
        const PetscInt gi_m = nodes->GI_U[m];

        Found = PETSC_FALSE;
        for(u=0; u<ln_S[n][1]; u++) 
          if(ln_S[n][u+2]==gi_m)  {
            Found = PETSC_TRUE;
            break;}
        if (!Found)
        { ln_S[n][1]++; ln_S[n][ln_S[n][1]+1] = gi_m; }
      }
    }
  }

  if (size > 1)
  {
    
    size_gn_M = 0, size_gn_S = 0;
    for (i=0; i < nodes->N_U_total_loc; i++) 
    {
      if (nodes->t[i] > 1) 
      { if (ln_M[i][1]>0) size_gn_M++;
        if (ln_S[i][1]>0) size_gn_S++;  }
    }
    
    ierr = PetscMalloc1(size_gn_M, &gn_M);
    for (i=0; i < size_gn_M; i++)
      ierr = PetscMalloc1(A_Width, &gn_M[i]), gn_M[i][0] = -1, gn_M[i][1] = 0; 
    ierr = PetscMalloc1(size_gn_S, &gn_S);
    for (i=0; i < size_gn_S; i++)
      ierr = PetscMalloc1(A_Width, &gn_S[i]), gn_S[i][0] = -1, gn_S[i][1] = 0; 

    k = 0, u = 0;
    for (i=0; i < nodes->N_U_total_loc; i++)
    {    
      if (nodes->t[i] < 2) 
        continue;   

      if  (ln_M[i][1]>0)
      {
        gn_M[k][0] = ln_M[i][0];            
        gn_M[k][1] = ln_M[i][1];            
        for (j = 0; j < ln_M[i][1]; j++)    
          gn_M[k][j+2] = ln_M[i][j+2];      
        k++;
      }

      if  (ln_S[i][1]>0)
      {
        gn_S[u][0] = ln_S[i][0];            
        gn_S[u][1] = ln_S[i][1];            
        for (j = 0; j < ln_S[i][1]; j++)    
          gn_S[u][j+2] = ln_S[i][j+2];      
        u++;
      }
    }
    
    for(r=rank; r>0; r--) 
    {
      ierr = MPI_Isend( &size_gn_M, 1, MPIU_INT, r-1, 0, PETSC_COMM_WORLD, &Request);         CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request);   CHKERRMPI(ierr);
      for(i=0; i<size_gn_M; i++)
      { ierr = MPI_Isend( gn_M[i], A_Width, MPIU_INT, r-1, i+1, PETSC_COMM_WORLD, &Request);  CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);}

      ierr = MPI_Isend( &size_gn_S, 1, MPIU_INT, r-1, 100000, PETSC_COMM_WORLD, &Request);  CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
      for(i=0; i<size_gn_S; i++)
      { ierr = MPI_Isend( gn_S[i], A_Width, MPIU_INT, r-1, i+100001, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request);     CHKERRMPI(ierr);}
    }
    
    for (r=rank+1; r<size; r++) 
    {
      
      flg = 0; 
      while(!flg) MPI_Iprobe(r, 0, PETSC_COMM_WORLD, &flg, &Status );
      ierr = MPI_Irecv( &rsize_gn_M, 1, MPIU_INT, r, 0, PETSC_COMM_WORLD, &Request);        CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
      ierr = PetscMalloc1(rsize_gn_M, &rgn_M); CHKERRQ(ierr);
      for (i=0; i < rsize_gn_M; i++){
        ierr = PetscMalloc1(A_Width, &rgn_M[i]); CHKERRQ(ierr);
        flg = 0; 
        while(!flg) MPI_Iprobe(r, i+1, PETSC_COMM_WORLD, &flg, &Status );
        ierr = MPI_Irecv(rgn_M[i], A_Width, MPIU_INT, r, i+1, PETSC_COMM_WORLD, &Request);   CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request);CHKERRMPI(ierr);}  

      flg = 0; 
      while(!flg) MPI_Iprobe(r, 100000, PETSC_COMM_WORLD, &flg, &Status );
      ierr = MPI_Irecv( &rsize_gn_S, 1, MPIU_INT, r, 100000, PETSC_COMM_WORLD, &Request);   CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
      ierr = PetscMalloc1(rsize_gn_S, &rgn_S); CHKERRQ(ierr);
      for (i=0; i < rsize_gn_S; i++){
        ierr = PetscMalloc1(A_Width, &rgn_S[i]); 
        flg = 0; 
        while(!flg) MPI_Iprobe(r, i+100001, PETSC_COMM_WORLD, &flg, &Status );
        ierr = MPI_Irecv(rgn_S[i], A_Width, MPIU_INT, r, i+100001, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);}
      
      for (i=0; i < rsize_gn_M; i++)
      {
        const PetscInt gi = rgn_M[i][0];
        if (gi<GI_start || gi >= GI_end)
          continue;
        else 
          for(li=0;li<nodes->N_U_total_loc;li++)
            if (nodes->GI_U[li] == gi)
              break;

        if(nodes->GI_U[li] != gi)
        { PetscErrorPrintf("Main global index %d was not found on rank %d.\n", gi, rank);
          SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Global index not found."); }

        for(j=0;j< rgn_M[i][1]; j++)
        {
          const PetscInt gi2 = rgn_M[i][2+j];

          Found = PETSC_FALSE;
          for(u=0; u<ln_M[li][1]; u++)
            if(ln_M[li][u+2]==gi2)  {
              Found = PETSC_TRUE;
              break;}
          if (!Found)
          { ln_M[li][1]++;
            
            ln_M[li][ln_M[li][1]+1] = gi2; }
        }
      }

      for (i=0; i < rsize_gn_S; i++)
      {
        const PetscInt gi = rgn_S[i][0];

        if (gi<GI_start || gi >= GI_end)
          continue;
        else 
          for(li=0;li<nodes->N_U_total_loc;li++)
            if (nodes->GI_U[li] == gi)
              break;

        if(nodes->GI_U[li] != gi)
        { PetscErrorPrintf("Side global index %d was not found on rank %d.\n", gi, rank);
          SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Global index not found."); }

        for(j=0;j< rgn_S[i][1]; j++)
        {
          const PetscInt gi2 = rgn_S[i][2+j];

          Found = PETSC_FALSE;
          for(u=0; u<ln_S[li][1]; u++)
            if(ln_S[li][u+2]==gi2)  {
              Found = PETSC_TRUE;
              break;}
          if (!Found)
          { ln_S[li][1]++; 
            
            ln_S[li][ln_S[li][1]+1] = gi2; }
        }
      }

      for (i=0; i < rsize_gn_M; i++) PetscFree(rgn_M[i]); 
      PetscFree(rgn_M); 
      for (i=0; i < rsize_gn_S; i++) PetscFree(rgn_S[i]); 
      PetscFree(rgn_S); 
    }

    for (i=0; i < size_gn_M; i++)    PetscFree(gn_M[i]);
    PetscFree(gn_M);
    for (i=0; i < size_gn_S; i++)    PetscFree(gn_S[i]);
    PetscFree(gn_S);
  }

  return 0;
}


PetscErrorCode CreateFullGlobalMatrix(Mat *A, const PetscInt glo, const PetscInt loc)
{
  PetscErrorCode ierr;
  ierr = MatCreate(PETSC_COMM_WORLD, A); CHKERRQ(ierr);
  ierr = MatSetSizes(*A, loc, loc, glo, glo); CHKERRQ(ierr);
  ierr = MatSetFromOptions(*A); CHKERRQ(ierr);
  ierr = MatSetUp(*A); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode CreateFullPreDeterminedStokesGlobalMatrix( Mat *A, System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements,
  StokesVariables_Struct *stokesVariables)
{  
  PetscErrorCode ierr;

  const PetscInt rank       = system->rank;
  const PetscInt Dim        = system->Dim;
  const PetscInt Nv_U       = stokesVariables->Nv_U;
  const PetscInt Nv_P       = stokesVariables->Nv_P;
  const PetscInt GI_V_start = stokesVariables->GI_V_start;
  const PetscInt glo        = stokesVariables->N_V_total_glo;
  const PetscInt loc        = stokesVariables->N_V_total_loc; 
  const PetscInt N_U_total_loc = nodes->N_U_total_loc;
  const PetscScalar zero = 0.0;

  PetscInt GI_start, GI_end; 

  PetscInt i, u, r;
  PetscInt *d_nnz, *o_nnz; 
  PetscInt *LI_V;          
  PetscInt Mn_l, Mn_n;     
  PetscInt Sn_l, Sn_n ;    

  PetscInt **ln_M,**ln_S;  

  PetscInt A_Width = 0;    
  
  if (Dim == 2)       A_Width = 30;
  else if (Dim == 3)  A_Width = 100;

  ierr = PetscMalloc2(N_U_total_loc, &ln_M, N_U_total_loc, &ln_S); CHKERRQ(ierr);
  for (i=0; i < N_U_total_loc; i++)
  { 
    PetscMalloc2(A_Width, &ln_M[i], A_Width, &ln_S[i]);
    ln_M[i][0] = -1, ln_M[i][1] = 0;
    ln_S[i][0] = -1, ln_S[i][1] = 0;  
  }

  CalculateTheNumberOfConnections(system, nodes, elements, ln_M, ln_S, A_Width);
  
  for (i=0; i<N_U_total_loc; i++)
  {
    if(ln_M[i][1]> A_Width-2)
    { PetscErrorPrintf("Error in rank %d for ln_M at row %d an Gi %d. Number of connections = %d > Max connection %d\n",rank, i, ln_M[i][0], ln_M[i][1], A_Width-2);
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! in ln_M"); }

    if(ln_S[i][1]> A_Width-2)
    { PetscErrorPrintf("Error in rank %d for ln_S at row %d an Gi %d. Number of connections = %d > Max connection %d\n",rank, i, ln_S[i][0], ln_M[i][1], A_Width-2);
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! in ln_S"); }
  }
 
  GI_start = 0;
  for (r=0; r < rank; r++)  GI_start += nodes->N_cores[r][2] - nodes->N_cores[r][3];
  GI_end = GI_start + nodes->N_cores[rank][2] - nodes->N_cores[rank][3];

  PetscMalloc1(Nv_U+Nv_P, &LI_V);
  PetscMalloc2(loc, &d_nnz, loc, &o_nnz);

  PetscPrintf(PETSC_COMM_WORLD,"  Building a predetermined global A matrix with space for penalty at cont. rows on all ranks.\n");

  for (i=0; i < N_U_total_loc; i++)
  {
    if (nodes->t[i] > 1) 
      continue;      
    
    for(u=0; u<Nv_U; u++)  LI_V[u] = stokesVariables->GI_V[u][i] - GI_V_start;

    Mn_l = 0, Mn_n = 0;
    for(u=0; u<ln_M[i][1]; u++)
    {
      const PetscInt gi = ln_M[i][u+2];
      if (gi<GI_start || gi >= GI_end) Mn_n++;
      else Mn_l++;
    }

    Sn_l = 0, Sn_n = 0;
    for(u=0; u<ln_S[i][1]; u++)
    {
      const PetscInt gi = ln_S[i][u+2];
      if (gi<GI_start || gi >= GI_end) Sn_n++;
      else Sn_l++;
    }

    if (nodes->LI_UtoP[i] != -1) 
    {
      
      for(u=Nv_U; u<Nv_U+Nv_P; u++) 
        LI_V[u] = stokesVariables->GI_V[u][nodes->LI_UtoP[i]] - GI_V_start;

      for(u=0; u<Nv_U; u++)
        d_nnz[LI_V[u]] = 2*Mn_l+Sn_l+2, 
        o_nnz[LI_V[u]] = 2*Mn_n+Sn_n+2;
      
      for(u=Nv_U; u<Nv_U+Nv_P; u++)
        d_nnz[LI_V[u]] = (Nv_U+Nv_P)*(Mn_l+Sn_l), 
        o_nnz[LI_V[u]] = (Nv_U+Nv_P)*(Mn_n+Sn_n);
    }
    else 
    {
      for(u=0; u<Nv_U; u++)
        d_nnz[LI_V[u]] = 2*Mn_l+Sn_l+2,
        o_nnz[LI_V[u]] = 2*Mn_n+Sn_n+2;
    }
  }

  for (i=0; i < loc; i++)
    if (o_nnz[i] < zero || d_nnz[i] < zero)
    { PetscErrorPrintf("o_nnz or d_nnz has a negative value of %d or %d at %d.\n",o_nnz[i], d_nnz[i], i);
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Incorrect mesh structure."); }

  PetscBarrier(NULL);

  for (i=0; i < N_U_total_loc; i++)
    ierr = PetscFree2(ln_M[i],ln_S[i]); CHKERRQ(ierr);
  ierr = PetscFree2(ln_M, ln_S); CHKERRQ(ierr);

  ierr = MatCreateAIJ(PETSC_COMM_WORLD, loc, loc, glo, glo, PETSC_DECIDE, d_nnz, PETSC_DECIDE, o_nnz, A);
  ierr = MatSetFromOptions(*A);
  ierr = MatSetUp(*A);

  ierr = PetscPrintf(PETSC_COMM_WORLD,"    Finished building a predetermined global matrix A of size %d on all ranks.\n", glo);

  PetscFree2(d_nnz, o_nnz);
  PetscFree(LI_V);
  return ierr;
}

#endif