
  typedef struct {

    PetscInt N_P_total_glo;
    PetscInt N_P_inter_glo;
    PetscInt N_P_bound_glo;
    PetscInt N_P_total_loc;
    PetscInt N_P_inter_loc;
    PetscInt N_P_bound_loc;
    PetscInt N_P_ghost_loc;
    PetscInt N_U_total_glo;
    PetscInt N_U_inter_glo;
    PetscInt N_U_bound_glo;
    PetscInt N_U_total_loc;
    PetscInt N_U_inter_loc;
    PetscInt N_U_bound_loc;
    PetscInt N_U_ghost_loc;
    PetscInt N_Gr;         
    PetscScalar *x, *y, *z; 
    PetscInt    *t;         
    PetscInt    **e;        
    PetscInt    *GI_U;          
    PetscInt    *LI_UtoP;       
    PetscInt    *LI_U_bound;    
    PetscInt    *Gr_U_bound;    
    PetscInt    *LI_U_ghost;    
    PetscInt    *Ra_U_ghost;    
    PetscInt    *Ra_LI_U_ghost; 
    PetscInt    *GI_P;          
    PetscInt    *LI_PtoU;       
    PetscInt    *LI_P_ghost;    
    PetscInt    *Ra_P_ghost;    
    PetscInt    *Ra_LI_P_ghost;  
    PetscInt    *Ra_ULI_P_ghost; 
    PetscInt   **N_cores;
    PetscInt     GI_U_start, GI_U_end;
    PetscInt     GI_P_start, GI_P_end;
    PetscScalar  DC[3]; 

  } Nodes_Struct; // Linear System


void NodesInitialize(Nodes_Struct *nodes)
{
  nodes->N_P_total_glo = -1;
  nodes->N_P_inter_glo = -1;
  nodes->N_P_bound_glo = -1;
  nodes->N_U_total_glo = -1;
  nodes->N_U_inter_glo = -1;
  nodes->N_U_bound_glo = -1;
  nodes->N_P_total_loc = -1;
  nodes->N_P_inter_loc = -1;
  nodes->N_P_bound_loc = -1;
  nodes->N_P_ghost_loc = -1;
  nodes->N_U_total_loc = -1;
  nodes->N_U_inter_loc = -1;
  nodes->N_U_bound_loc = -1;
  nodes->N_U_ghost_loc = -1;
  nodes->N_Gr = -1;  
}

PetscErrorCode NodesCountU(Nodes_Struct *nodes, const PetscMPIInt rank)
{
  nodes->N_U_inter_loc = 0;
  nodes->N_U_bound_loc = 0;
  nodes->N_U_ghost_loc = 0;

  if(nodes->N_U_total_loc <= 0)  {
    PetscErrorPrintf("Error N_U_total_loc = %d in Nodes_Struct object rank %d can not be <= 0\n", nodes->N_U_total_loc, rank);
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error nodes->N_U_total_loc");}

  for(PetscInt i=0; i<nodes->N_U_total_loc; ++i)  {
    if      (nodes->t[i]==0)                    nodes->N_U_inter_loc++;
    else if (nodes->t[i]==1)                    nodes->N_U_bound_loc++;
    else if (nodes->t[i]==2 || nodes->t[i]==3)  nodes->N_U_ghost_loc++;}

  if(nodes->N_U_total_loc != (nodes->N_U_inter_loc+nodes->N_U_bound_loc+nodes->N_U_ghost_loc))  {
    PetscErrorPrintf("Error N_U_total_loc = %d different from (N_U_inter_loc+N_U_bound_loc+N_U_ghost_loc) = %d in Nodes_Struct object on rank %d\n", 
        nodes->N_U_total_loc, (nodes->N_U_inter_loc+nodes->N_U_bound_loc+nodes->N_U_ghost_loc), rank);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error nodes->N_U_total_loc"); }
  return 0; 
}
PetscErrorCode NodesMallocUArraysForAll(Nodes_Struct *nodes, const PetscMPIInt size, const PetscMPIInt rank) 
{
  if(nodes->N_U_total_loc <= 0)  {
    PetscErrorPrintf("Error N_U_total_loc = %d in Nodes_Struct object rank %d can not be <= 0\n", nodes->N_U_total_loc, rank);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error nodes->N_U_total_loc");      }
  
  PetscMalloc6(nodes->N_U_total_loc, &nodes->x, nodes->N_U_total_loc, &nodes->y, nodes->N_U_total_loc, &nodes->z, 
               nodes->N_U_total_loc, &nodes->t, nodes->N_U_total_loc, &nodes->GI_U, nodes->N_U_total_loc, &nodes->LI_UtoP);
  PetscMalloc1(nodes->N_U_total_loc, &nodes->e);
  for(PetscInt i=0; i<nodes->N_U_total_loc; i++) {
    nodes->GI_U[i] = -1;
    PetscMalloc1(3, &nodes->e[i]);
    nodes->e[i][0]=0, nodes->e[i][1]=0, nodes->e[i][2]=0, nodes->LI_UtoP[i] = -1;}
  PetscMalloc1(size, &nodes->N_cores);
  for(PetscInt i=0; i<size; i++)  PetscMalloc1(4, &nodes->N_cores[i]);
  return 0; 
}

PetscErrorCode NodesMallocPArraysForAll(Nodes_Struct *nodes, const PetscMPIInt rank) 
{
  if(nodes->N_P_total_loc <= 0)  {
    PetscErrorPrintf("Error N_P_total_loc  = %d in Nodes_Struct object rank %d can not be <= 0\n", nodes->N_P_total_loc, rank);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error nodes->N_P_total_loc");      }
  PetscMalloc2(nodes->N_P_total_loc, &nodes->GI_P, nodes->N_P_total_loc, &nodes->LI_PtoU);
  return 0; 
}

PetscErrorCode NodesMallocUBoundaryArrays(Nodes_Struct *nodes, const PetscMPIInt rank)  
{
  if(nodes->N_U_bound_loc < 0)  {
    PetscErrorPrintf("Error N_U_bound_loc = %d in Nodes_Struct object rank %d can not be <= 0\n", nodes->N_U_bound_loc, rank);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error nodes->N_U_bound_loc");}      
  PetscMalloc2(nodes->N_U_bound_loc, &nodes->LI_U_bound, nodes->N_U_bound_loc, &nodes->Gr_U_bound);
  return 0; 
}

PetscErrorCode NodesMallocUGhostArrays(Nodes_Struct *nodes, const PetscMPIInt rank)  
{
  PetscMalloc3(nodes->N_U_ghost_loc, &nodes->LI_U_ghost, nodes->N_U_ghost_loc, &nodes->Ra_U_ghost, nodes->N_U_ghost_loc, &nodes->Ra_LI_U_ghost);
  for (PetscInt i=0; i < nodes->N_U_ghost_loc; i++)
    nodes->LI_U_ghost[i] = 0, nodes->Ra_U_ghost[i] = 0, nodes->Ra_LI_U_ghost[i] = 0;
  return 0; 
}

PetscErrorCode NodesMallocPGhostArrays(Nodes_Struct *nodes, const PetscMPIInt rank)  
{
  PetscMalloc4(nodes->N_P_ghost_loc, &nodes->LI_P_ghost, nodes->N_P_ghost_loc, &nodes->Ra_P_ghost, 
    nodes->N_P_ghost_loc, &nodes->Ra_LI_P_ghost, nodes->N_P_ghost_loc, &nodes->Ra_ULI_P_ghost);
  for (PetscInt i=0; i < nodes->N_P_ghost_loc; i++)
    nodes->LI_P_ghost[i] = 0, nodes->Ra_P_ghost[i] = 0, nodes->Ra_LI_P_ghost[i] = 0, nodes->Ra_ULI_P_ghost[i] = 0;
  return 0; 
}

PetscErrorCode NodesFree(Nodes_Struct *nodes, const PetscMPIInt size) 
{
  PetscPrintf(PETSC_COMM_WORLD,"  Nodes memory is being freed\n");
  for(PetscInt i=0; i<nodes->N_U_total_loc; i++) PetscFree(nodes->e[i]);
  PetscFree(nodes->e);
  PetscFree6(nodes->x, nodes->y, nodes->z, nodes->t, nodes->GI_U, nodes->LI_UtoP);
  for(PetscInt i=0; i<size; i++)  PetscFree(nodes->N_cores[i]);
  PetscFree(nodes->N_cores);
  PetscFree2(nodes->GI_P, nodes->LI_PtoU);
  PetscFree2(nodes->LI_U_bound, nodes->Gr_U_bound);
  PetscFree3(nodes->LI_U_ghost, nodes->Ra_U_ghost, nodes->Ra_LI_U_ghost);
  PetscFree4(nodes->LI_P_ghost, nodes->Ra_P_ghost, nodes->Ra_LI_P_ghost, nodes->Ra_ULI_P_ghost);

  return 0; 
}

PetscErrorCode NodesPrintInfoMini(System_Struct *system, Nodes_Struct *nodes, PetscBool boolNodes) 
{
  if (boolNodes)
  {
    PetscPrintf(PETSC_COMM_WORLD,"  Nodes Info:\n");
    PetscPrintf(PETSC_COMM_WORLD,"  Total global U nodes: %d\n", nodes->N_U_total_glo);
    PetscPrintf(PETSC_COMM_WORLD,"  Total global internal U nodes: %d\n", nodes->N_U_inter_glo);
    PetscPrintf(PETSC_COMM_WORLD,"  Total global boundary U nodes: %d\n", nodes->N_U_bound_glo);

    PetscPrintf(PETSC_COMM_WORLD,"  Total local U nodes:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: %d\n", system->rank, nodes->N_U_total_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local internal U nodes:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: %d\n", system->rank, nodes->N_U_inter_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local boundary U nodes:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: %d\n", system->rank, nodes->N_U_bound_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local ghost U nodes:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: %d\n", system->rank, nodes->N_U_ghost_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

    PetscPrintf(PETSC_COMM_WORLD,"  Total global P nodes: %d\n", nodes->N_P_total_glo);
    PetscPrintf(PETSC_COMM_WORLD,"  Total global internal P nodes: %d\n", nodes->N_P_inter_glo);
    PetscPrintf(PETSC_COMM_WORLD,"  Total global boundary P nodes: %d\n", nodes->N_P_bound_glo);

    PetscPrintf(PETSC_COMM_WORLD,"  Total local P nodes:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: %d\n", system->rank, nodes->N_P_total_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local internal P nodes:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: %d\n", system->rank, nodes->N_P_inter_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local boundary P nodes:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: %d\n", system->rank, nodes->N_P_bound_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local ghost P nodes:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: %d\n", system->rank, nodes->N_P_ghost_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }
  return 0;
}
PetscErrorCode NodesPrintInfo(System_Struct *system, Nodes_Struct *nodes, PetscBool boolNodes) 
{  
  if (boolNodes)
  {
    NodesPrintInfoMini(system, nodes, boolNodes);

    PetscPrintf(PETSC_COMM_WORLD,"  x:\n");
    PetscScalarView(nodes->N_U_total_loc, nodes->x, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"  y:\n");
    PetscScalarView(nodes->N_U_total_loc, nodes->y, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"  z:\n");
    PetscScalarView(nodes->N_U_total_loc, nodes->z, stdout_viewer);

    PetscPrintf(PETSC_COMM_WORLD,"  t:\n");
    PetscIntView(nodes->N_U_total_loc, nodes->t, stdout_viewer);

    PetscPrintf(PETSC_COMM_WORLD,"  GI_U:\n");
    PetscIntView(nodes->N_U_total_loc, nodes->GI_U, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"  LI_UtoP:\n");
    PetscIntView(nodes->N_U_total_loc, nodes->LI_UtoP, stdout_viewer);

    PetscPrintf(PETSC_COMM_WORLD,"  LI_U_bound:\n");
    PetscIntView(nodes->N_U_bound_loc, nodes->LI_U_bound, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"  Gr_U_bound:\n");
    PetscIntView(nodes->N_U_bound_loc, nodes->Gr_U_bound, stdout_viewer);

    PetscPrintf(PETSC_COMM_WORLD,"  LI_U_ghost:\n");
    PetscIntView(nodes->N_U_ghost_loc, nodes->LI_U_ghost, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"  Ra_U_ghost:\n");
    PetscIntView(nodes->N_U_ghost_loc, nodes->Ra_U_ghost, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"  Ra_LI_U_ghost:\n");
    PetscIntView(nodes->N_U_ghost_loc, nodes->Ra_LI_U_ghost, stdout_viewer);

    PetscPrintf(PETSC_COMM_WORLD,"  GI_P:\n");
    PetscIntView(nodes->N_P_total_loc, nodes->GI_P, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"  LI_UtoP:\n");
    PetscIntView(nodes->N_P_total_loc, nodes->LI_UtoP, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"  LI_P_ghost:\n");
    PetscIntView(nodes->N_P_ghost_loc, nodes->LI_P_ghost, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"  Ra_P_ghost:\n");
    PetscIntView(nodes->N_P_ghost_loc, nodes->Ra_P_ghost, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"  Ra_LI_P_ghost:\n");
    PetscIntView(nodes->N_P_ghost_loc, nodes->Ra_LI_P_ghost, stdout_viewer);
  }
  return 0; 
}

PetscInt NodesFindTheMaxGI(System_Struct *system, Nodes_Struct *nodes)
{
  PetscInt local_max_GI;       
  PetscInt global_max_GI;    
  PetscMPIInt size = system->size;
  PetscMPIInt rank = system->rank;     
  local_max_GI= -1000000;
  for(PetscInt i=0; i<nodes->N_U_total_loc; i++)
    if (nodes->GI_U[i] > local_max_GI)
      local_max_GI = nodes->GI_U[i];

  if (size == 1)      global_max_GI = local_max_GI;
  else                MaxIntAcrossRanks(size, rank, &global_max_GI, local_max_GI);

  return global_max_GI;
}



PetscErrorCode NodesSearch(System_Struct *system, Nodes_Struct *nodes, PetscScalar *NS, 
  PetscInt *Li, PetscInt *rankFound)
{
  PetscErrorCode ierr;
  PetscScalar Dmin = 1e8, D;
  PetscInt li;
  MPI_Request Request;
  MPI_Status  Status;

  *Li=-1;
  *rankFound=-1;

  for (li=0; li<nodes->N_U_total_loc; li++)
  {
    D = CalcDistance2(nodes->x[li], nodes->y[li], nodes->z[li], NS[0], NS[1], NS[2]);
    if (D<Dmin && nodes->t[li]==1)
    {
      Dmin = D;
      *Li = li;
      *rankFound = system->rank;
    }
  } 

  // Compare the Ds of all ranks to see who is the smallest
  if (system->size>1)
  { 
    if (system->rank != 0) // Send to 0 if found or not
    {
      ierr = MPI_Isend(&Dmin, 1, MPIU_SCALAR, 0, system->rank+10000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
    }
    else // Rank == 0
    {
      for(PetscInt r = 1; r<system->size; ++r) // r is the rank that will send points
      {
        ierr = MPI_Irecv(&D, 1, MPIU_SCALAR,  r, r+10000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
        if (D< Dmin)
          Dmin = D, *rankFound = r;
      }
    }
    ierr = MPI_Bcast( rankFound, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); // Send the GI_P row that will be zerod to all
  }

  if (*rankFound==-1){
    PetscErrorPrintf("Error! The node was not found \n");
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in  NodesSearch"); }
  else
    PetscPrintf(PETSC_COMM_WORLD, "Node was found on rank %d\n", *rankFound);

  return 0;
}
