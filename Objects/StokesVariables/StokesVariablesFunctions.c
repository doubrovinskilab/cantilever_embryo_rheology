#ifndef StokesVariables_Functions_C
#define StokesVariables_Functions_C

void StokesVariablesInitialize(StokesVariables_Struct *stokesVariables)
{
  stokesVariables->Nv = -1;
  stokesVariables->Nv_U = -1;
  stokesVariables->Nv_P = -1;
  stokesVariables->N_Gr = -1;
  stokesVariables->N_V_total_glo = -1;
  stokesVariables->N_V_total_loc = -1;
  stokesVariables->N_V_inter_loc = -1;
  stokesVariables->N_V_bound_loc = -1;
  stokesVariables->N_V_ghost_loc = -1;
  stokesVariables->N_V_diric_loc = -1;
  stokesVariables->N_V_neuma_loc = -1;
  stokesVariables->N_V_robin_loc = -1;
  stokesVariables->N_V_unkno_loc = -1;
  stokesVariables->N_V_U_total_glo = -1;
  stokesVariables->N_V_P_total_glo = -1;
  stokesVariables->N_V_U_total_loc = -1;
  stokesVariables->N_V_P_total_loc = -1;
  stokesVariables->N_V_U_diric_loc = -1;
  stokesVariables->N_V_P_diric_loc = -1;
  stokesVariables->GI_V_start = -1;
  stokesVariables->GI_V_P_start = -1;
  stokesVariables->N_Gr = -1;
  stokesVariables->ViscMethod = ONE_CONSTANT;
  stokesVariables->mu2 = -1;
  stokesVariables->A = -1;
  stokesVariables->B = -1;
  stokesVariables->C = -1;
  stokesVariables->Bd = -1;
  stokesVariables->Bv = -1;
}


PetscErrorCode StokesVariablesFree(StokesVariables_Struct *stokesVariables) 
{
  PetscErrorCode ierr;
  for(PetscInt i = 0; i < stokesVariables->Nv; ++i){
    ierr = PetscFree(stokesVariables->GI_V[i]); CHKERRQ(ierr);
    ierr = PetscFree3(stokesVariables->Gr_bound_value[i],stokesVariables->Gr_bound_value2[i],stokesVariables->Gr_bound_type[i]); CHKERRQ(ierr);}
  ierr = PetscFree(stokesVariables->GI_V); CHKERRQ(ierr);
  ierr = PetscFree3(stokesVariables->Gr_bound_value, stokesVariables->Gr_bound_value2, stokesVariables->Gr_bound_type); CHKERRQ(ierr);
  ierr = PetscFree(stokesVariables->GI_V_ghost); CHKERRQ(ierr);
  ierr = PetscFree(stokesVariables->GI_V_U_bound_dir); CHKERRQ(ierr);
  ierr = PetscFree(stokesVariables->GI_V_P_bound_dir); CHKERRQ(ierr);
  ierr = PetscFree(stokesVariables->GI_V_U_non_dir); CHKERRQ(ierr);
  ierr = PetscFree(stokesVariables->GI_V_P_non_dir); CHKERRQ(ierr);
  return 0; 
}

PetscErrorCode StokesVariablesBuild(System_Struct *system, Nodes_Struct *nodes, StokesVariables_Struct *stokesVariables) 
{
  PetscErrorCode ierr;
  PetscInt i, j;
  stokesVariables->N_Gr = nodes->N_Gr;
  if(system->Dim == 2){
    stokesVariables->Nv   = 3;
    stokesVariables->Nv_U = 2;
    stokesVariables->Nv_P = 1;}
  else if(system->Dim == 3){
    stokesVariables->Nv   = 4;
    stokesVariables->Nv_U = 3;
    stokesVariables->Nv_P = 1;}
   else {
    PetscErrorPrintf("Error in Variables_Class. Dimenstion should be 2 or 3.\n");
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong dimension size."); }

  const PetscInt N_Gr = stokesVariables->N_Gr;
  const PetscInt Nv   = stokesVariables->Nv;
  const PetscInt Nv_P = stokesVariables->Nv_P;
  const PetscInt Nv_U = stokesVariables->Nv_U;

  stokesVariables->N_P_total_loc = nodes->N_P_total_loc;
  stokesVariables->N_U_total_loc = nodes->N_U_total_loc;

  stokesVariables->N_V_total_glo = Nv_P*nodes->N_P_total_glo + Nv_U*nodes->N_U_total_glo; 

  stokesVariables->N_V_ghost_loc = Nv_P*nodes->N_P_ghost_loc + Nv_U*nodes->N_U_ghost_loc;
  stokesVariables->N_V_bound_loc = Nv_P*nodes->N_P_bound_loc + Nv_U*nodes->N_U_bound_loc;
  stokesVariables->N_V_inter_loc = Nv_P*nodes->N_P_inter_loc + Nv_U*nodes->N_U_inter_loc;
  stokesVariables->N_V_total_loc = Nv_P*nodes->N_P_total_loc + Nv_U*nodes->N_U_total_loc - stokesVariables->N_V_ghost_loc;

  stokesVariables->N_V_P_total_loc = Nv_P*(nodes->N_P_total_loc - nodes->N_P_ghost_loc);
  stokesVariables->N_V_U_total_loc = Nv_U*(nodes->N_U_total_loc - nodes->N_U_ghost_loc);

  stokesVariables->N_V_U_total_glo = Nv_U*nodes->N_U_total_glo;
  stokesVariables->N_V_P_total_glo = Nv_P*nodes->N_P_total_glo;

  stokesVariables->GI_V_start = 0;
  for (i=0; i<system->rank; i++) 
    stokesVariables->GI_V_start += (nodes->N_cores[i][2]- nodes->N_cores[i][3])*Nv_U + (nodes->N_cores[i][0]-nodes->N_cores[i][1])*Nv_P;
  stokesVariables->GI_V_P_start = stokesVariables->GI_V_start + (nodes->N_cores[i][2]- nodes->N_cores[i][3])*Nv_U;

  ierr = PetscMalloc3(Nv, &stokesVariables->Gr_bound_value, Nv, &stokesVariables->Gr_bound_value2, Nv, &stokesVariables->Gr_bound_type); CHKERRQ(ierr);
  for(i = 0; i < Nv; ++i) 
  { ierr = PetscMalloc3(N_Gr, &stokesVariables->Gr_bound_value[i], N_Gr, &stokesVariables->Gr_bound_value2[i], N_Gr, &stokesVariables->Gr_bound_type[i]);
    CHKERRQ(ierr); }

  ierr = PetscMalloc1(Nv, &stokesVariables->GI_V); 

  for(i = 0; i < Nv_U; ++i) {
    ierr = PetscMalloc1(stokesVariables->N_U_total_loc, &stokesVariables->GI_V[i]);   CHKERRQ(ierr); 
    for(j=0; j < stokesVariables->N_U_total_loc; j++) stokesVariables->GI_V[i][j] = -1;}
  
  ierr = PetscMalloc1(stokesVariables->N_P_total_loc, &stokesVariables->GI_V[Nv-1]);  CHKERRQ(ierr); 
  for(j=0; j < stokesVariables->N_P_total_loc; j++) stokesVariables->GI_V[Nv-1][j] = -1;

  ierr = PetscMalloc1(stokesVariables->N_V_ghost_loc, &stokesVariables->GI_V_ghost); CHKERRQ(ierr);
  ierr = PetscBarrier(NULL); CHKERRQ(ierr);
  PetscInt GI_V_start_r, GI_V_Delta_r, GI_Largest;
  PetscInt v, r;
  GI_V_start_r = stokesVariables->GI_V_start;
  i = 0;
  for(v =0; v < Nv_U; ++v) 
    for (PetscInt n =0; n < nodes->N_U_total_loc; ++n)
      if (nodes->t[n] < 2)          
        stokesVariables->GI_V[v][n] = GI_V_start_r + i, i++;
      else                          
        stokesVariables->GI_V[v][n] = -1;
  for (PetscInt n =0; n < nodes->N_P_total_loc; ++n)
    if (nodes->t[nodes->LI_PtoU[n]] < 2)    
      stokesVariables->GI_V[Nv_U][n] = GI_V_start_r + i, i++;
    else                                    
      stokesVariables->GI_V[Nv_U][n] = -1;

  j = 0;
  for(i=0; i < nodes->N_U_ghost_loc; i++)
  {
    const PetscInt Li   = nodes->LI_U_ghost[i];
    const PetscInt GI_U = nodes->GI_U[Li];

    GI_Largest = 0, GI_V_start_r=0, GI_V_Delta_r=0;
    
    for(r=0; r<system->size; r++)
    {
      GI_Largest += nodes->N_cores[r][2] - nodes->N_cores[r][3];
      if (GI_U < GI_Largest)
        break; 
      GI_V_start_r += (nodes->N_cores[r][2]- nodes->N_cores[r][3])*Nv_U + (nodes->N_cores[r][0]-nodes->N_cores[r][1])*Nv_P;
    }
    GI_V_Delta_r = GI_U - (GI_Largest - (nodes->N_cores[r][2] - nodes->N_cores[r][3]));

    stokesVariables->GI_V[0][Li] = GI_V_start_r + GI_V_Delta_r;
    stokesVariables->GI_V[1][Li] = GI_V_start_r + GI_V_Delta_r + nodes->N_cores[r][2] - nodes->N_cores[r][3] ;
    if (system->Dim == 3)
      stokesVariables->GI_V[2][Li] = GI_V_start_r + GI_V_Delta_r + 2*(nodes->N_cores[r][2]-nodes->N_cores[r][3]);

    
    stokesVariables->GI_V_ghost[j]  = stokesVariables->GI_V[0][Li];
    stokesVariables->GI_V_ghost[j+nodes->N_U_ghost_loc] = stokesVariables->GI_V[1][Li];
    if (system->Dim == 3) 
      stokesVariables->GI_V_ghost[j+2*nodes->N_U_ghost_loc] = stokesVariables->GI_V[2][Li];
    j++;
  }

  
  j = 0;
  for(i=0; i < nodes->N_P_ghost_loc; i++)
  {
    const PetscInt Li   = nodes->LI_P_ghost[i];
    const PetscInt GI_P = nodes->GI_P[Li];

    GI_Largest = 0, GI_V_start_r=0, GI_V_Delta_r=0;
    
    for(r=0; r<system->size; r++)
    {
      GI_Largest += nodes->N_cores[r][0] - nodes->N_cores[r][1];
      if (GI_P < GI_Largest)
        break; 
      GI_V_start_r += (nodes->N_cores[r][2]- nodes->N_cores[r][3])*Nv_U + (nodes->N_cores[r][0]-nodes->N_cores[r][1])*Nv_P;
    } 
    GI_V_start_r += (nodes->N_cores[r][2]- nodes->N_cores[r][3])*Nv_U;
    GI_V_Delta_r  = GI_P - (GI_Largest - (nodes->N_cores[r][0] - nodes->N_cores[r][1]));
    
    if (system->Dim == 2)
      stokesVariables->GI_V[2][Li] = GI_V_start_r + GI_V_Delta_r, 
      stokesVariables->GI_V_ghost[j+2*nodes->N_U_ghost_loc] = stokesVariables->GI_V[2][Li];
    else if (system->Dim == 3)
      stokesVariables->GI_V[3][Li] = GI_V_start_r + GI_V_Delta_r, 
      stokesVariables->GI_V_ghost[j+3*nodes->N_U_ghost_loc] = stokesVariables->GI_V[3][Li];
    j++;
  }

  ierr = PetscBarrier(NULL); CHKERRQ(ierr);
  return 0; 
}

PetscErrorCode StokesVariablesReadBoundary(System_Struct *system, Nodes_Struct *nodes, StokesVariables_Struct *stokesVariables) 
{
  PetscErrorCode ierr;
  char **data;
  char str_boundary[256]="";
  char str_value[256]="";
  
  PetscBool flg_i, flg_d, flg_n, flg_r, flg_s, flg_se;
  PetscInt N = 0;
  PetscInt n, v, j, g;
  PetscScalar value;

  for (v=0; v<stokesVariables->Nv; v++)
  {
    
    if (v == 0)
      PetscStrcpy(str_boundary,"Boundary_Ux"); 
    else if (v==1)
      PetscStrcpy(str_boundary,"Boundary_Uy"); 
    else if (v==2 && system->Dim == 3)
      PetscStrcpy(str_boundary,"Boundary_Uz"); 
    else if ((v==2 && system->Dim == 2) || (v==3 && system->Dim == 3))
      PetscStrcpy(str_boundary,"Boundary_P"); 
      
    N = SearchInputAsciiFileForNumberOfOccurances(system->SYSInputFile, str_boundary);
    if (N != stokesVariables->N_Gr)  {
      PetscErrorPrintf("The number of groups for fluid in InputFile for '%s' is %d is different from the number of groups in the mesh file which is %i\n",
      str_boundary,  N, stokesVariables->N_Gr);
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! number of groups in InputFile");  }

    for (n=0; n<N; n++)
    {
      SearchInputAsciiFileForNthOccurance(system->SYSInputFile, str_boundary, str_value, n);
      PetscStrToArray(str_value,',',&j, &data);

      ierr = PetscOptionsStringToInt(data[0], &g); CHKERRQ(ierr);
      PetscStrcmp(data[1], "i", &flg_i); 
      PetscStrcmp(data[1], "d", &flg_d); 
      PetscStrcmp(data[1], "n", &flg_n); 
      PetscStrcmp(data[1], "r", &flg_r); 
      PetscStrcmp(data[1], "s", &flg_s); 
      PetscStrcmp(data[1], "es",&flg_se);
      if(flg_i) 
      { 
        stokesVariables->Gr_bound_type[v][g]   = INTERNAL;
        if (j != 2)
          SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong fluid boundary."); 
        stokesVariables->Gr_bound_value[v][g]  = -1;  
        stokesVariables->Gr_bound_value2[v][g] = -1;  
        PetscPrintf(PETSC_COMM_WORLD,"  Fluid %s group %d: %s\n", str_boundary, g, FLUID_BC_TYPE_CHAR[stokesVariables->Gr_bound_type[v][g]] );
      }
      else if(flg_d) 
      { 
        stokesVariables->Gr_bound_type[v][g]  = DIRICHLET;
        if (j != 3)
          SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong fluid boundary."); 
        ierr = PetscOptionsStringToScalar(data[2], &value); CHKERRQ(ierr);
        stokesVariables->Gr_bound_value[v][g] = value;
        stokesVariables->Gr_bound_value2[v][g] = -1;  
        PetscPrintf(PETSC_COMM_WORLD,"  Fluid %s group %d: %s\n    v=%g\n", str_boundary, g, FLUID_BC_TYPE_CHAR[stokesVariables->Gr_bound_type[v][g]], stokesVariables->Gr_bound_value[v][g] );
      }
      else if(flg_n) 
      { 
        stokesVariables->Gr_bound_type[v][g]  = NEUMANN;
        if (j != 3) 
          SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong fluid boundary."); 
        ierr = PetscOptionsStringToScalar(data[2], &value); CHKERRQ(ierr);
        stokesVariables->Gr_bound_value[v][g] = value;
        stokesVariables->Gr_bound_value2[v][g] = -1;  
        PetscPrintf(PETSC_COMM_WORLD,"  Fluid %s group %d: %s\n    v'=B with B=%g\n", str_boundary, g, FLUID_BC_TYPE_CHAR[stokesVariables->Gr_bound_type[v][g]], stokesVariables->Gr_bound_value[v][g]);
      }      
      else if(flg_r) 
      { 
        stokesVariables->Gr_bound_type[v][g]  = ROBIN;
        if (j != 4)
          SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong fluid boundary."); 
        ierr = PetscOptionsStringToScalar(data[2], &value); CHKERRQ(ierr);
        stokesVariables->Gr_bound_value[v][g] = value;
        ierr = PetscOptionsStringToScalar(data[3], &value); CHKERRQ(ierr);
        stokesVariables->Gr_bound_value2[v][g] = value; 
        PetscPrintf(PETSC_COMM_WORLD,"  Fluid %s group %d: %s\n    v'+A*v=B with A=%g & B=%g\n", str_boundary, g, FLUID_BC_TYPE_CHAR[stokesVariables->Gr_bound_type[v][g]], stokesVariables->Gr_bound_value2[v][g],stokesVariables->Gr_bound_value[v][g]);
      }
      else if(flg_s) 
      { 
        stokesVariables->Gr_bound_type[v][g]  = SLIP;
        if (j != 4)
          SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong fluid boundary."); 
        ierr = PetscOptionsStringToScalar(data[2], &value); CHKERRQ(ierr);
        stokesVariables->Gr_bound_value[v][g] = value;
        ierr = PetscOptionsStringToScalar(data[3], &value); CHKERRQ(ierr);
        stokesVariables->Gr_bound_value2[v][g] = value; 
        PetscPrintf(PETSC_COMM_WORLD,"  Fluid %s group %d: %s\n    dv/dt+A*v=B with A=%g & B=%g and vn = 0\n", str_boundary, g, FLUID_BC_TYPE_CHAR[stokesVariables->Gr_bound_type[v][g]], stokesVariables->Gr_bound_value2[v][g],stokesVariables->Gr_bound_value[v][g]);
      }
      else if (flg_se)
      {
        stokesVariables->Gr_bound_type[v][g]  = EXP_SLIP;
        if (j != 4)
          SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong fluid boundary."); 
        ierr = PetscOptionsStringToScalar(data[2], &value); CHKERRQ(ierr);
        stokesVariables->Gr_bound_value[v][g] = value;
        ierr = PetscOptionsStringToScalar(data[3], &value); CHKERRQ(ierr);
        stokesVariables->Gr_bound_value2[v][g] = value; 
        PetscPrintf(PETSC_COMM_WORLD,"  Fluid %s group %d: %s\n    dv/dt+A*v=B with A=%g & B=%g and vn = 0\n", str_boundary, g, FLUID_BC_TYPE_CHAR[stokesVariables->Gr_bound_type[v][g]], stokesVariables->Gr_bound_value2[v][g],stokesVariables->Gr_bound_value[v][g]);        
      }
      else
        SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong fluid boundary."); 

      PetscStrToArrayDestroy(j, data);
    }
  }
  return 0; 
}


PetscErrorCode StokesVariableseSearchForDirichlet(Nodes_Struct *nodes, Elements_Struct *elements, StokesVariables_Struct *stokesVariables)
{
  PetscErrorCode  ierr;

  PetscInt i, j, k, v, q, gi_dir;

  const PetscInt Nv_U = stokesVariables->Nv_U;

  stokesVariables->N_V_diric_loc = 0;

  PetscInt Nb = nodes->N_U_bound_loc+ nodes->N_U_ghost_loc;

  PetscBool Found;
  
  PetscInt *temp_GI_V_U_bound_dir;
  ierr = PetscMalloc1(Nb*Nv_U, &temp_GI_V_U_bound_dir); CHKERRQ(ierr);
  for(i=0; i<Nb*Nv_U; i++) temp_GI_V_U_bound_dir[i] = -1;

  PetscInt *temp_GI_V_P_bound_dir;
  ierr = PetscMalloc1(Nb, &temp_GI_V_P_bound_dir); CHKERRQ(ierr);
  for(i=0; i<Nb; i++) temp_GI_V_P_bound_dir[i] = -1;
 
  j = 0, k = 0;
  for (PetscInt be=0; be < elements->N_bound_loc; ++be)
  {
    const PetscInt gr = elements->Gr_bound[be];

    
    for(v=0; v < Nv_U; v++)
    {
      if (stokesVariables->Gr_bound_type[v][gr] == DIRICHLET){ 
        for (i = 0; i < elements->U_S_bound; i++){
          const PetscInt li = elements->bound[be][i];
          const PetscInt gi = stokesVariables->GI_V[v][li];
          Found = PETSC_FALSE;
          for ( q = 0; q < j ; q++){
            if ( temp_GI_V_U_bound_dir[q] == gi){
              Found = PETSC_TRUE;
              break;
            }
          }
          if (PetscNot(Found) && gi>= stokesVariables->GI_V_start)  temp_GI_V_U_bound_dir[j] = gi, j++;
        }
      }
    }
    
    if (stokesVariables->Gr_bound_type[Nv_U][gr] == DIRICHLET) 
      for (i = 0; i < elements->P_S_bound; i++)
      {
        const PetscInt li_p = nodes->LI_UtoP[elements->bound[be][i]];
        const PetscInt gi_p = stokesVariables->GI_V[Nv_U][li_p];
        Found = PETSC_FALSE;
        for ( q = 0; q < k ; q++){
          if ( temp_GI_V_P_bound_dir[q] == gi_p){
            Found = PETSC_TRUE;
            break;
          }
        }
        if (PetscNot(Found) && gi_p>= stokesVariables->GI_V_P_start)  temp_GI_V_P_bound_dir[k] = gi_p, k++;
      }
  }

  if (j > Nb*Nv_U) {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong SearchForStokesDirichletVariables"); }

  if (k > Nb) {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong SearchForStokesDirichletVariables"); }

  quickSort(temp_GI_V_U_bound_dir, 0, j - 1);
  quickSort(temp_GI_V_P_bound_dir, 0, k - 1);

  stokesVariables->N_V_U_diric_loc = j; 
  ierr = PetscMalloc1(stokesVariables->N_V_U_diric_loc, &stokesVariables->GI_V_U_bound_dir); CHKERRQ(ierr);

  stokesVariables->N_V_P_diric_loc = k; 
  ierr = PetscMalloc1(stokesVariables->N_V_P_diric_loc, &stokesVariables->GI_V_P_bound_dir); CHKERRQ(ierr);

  stokesVariables->N_V_diric_loc = j + k; 
  stokesVariables->N_V_unkno_loc = stokesVariables->N_V_total_loc - stokesVariables->N_V_diric_loc;

  if (stokesVariables->N_V_diric_loc > stokesVariables->N_V_bound_loc){
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong SearchForStokesDirichletVariables"); }

  for (i=0; i < stokesVariables->N_V_U_diric_loc; i++)
    stokesVariables->GI_V_U_bound_dir[i] = temp_GI_V_U_bound_dir[i];

  for (i=0; i < stokesVariables->N_V_P_diric_loc; i++)
    stokesVariables->GI_V_P_bound_dir[i] = temp_GI_V_P_bound_dir[i];

  
  ierr = PetscFree(temp_GI_V_U_bound_dir); CHKERRQ(ierr);
  ierr = PetscFree(temp_GI_V_P_bound_dir); CHKERRQ(ierr);
  
  ierr = PetscMalloc1(stokesVariables->N_V_U_total_loc - stokesVariables->N_V_U_diric_loc, &stokesVariables->GI_V_U_non_dir); CHKERRQ(ierr);
  if (stokesVariables->N_V_U_diric_loc != 0)
  {
    i = 0, j = 0, gi_dir = stokesVariables->GI_V_U_bound_dir[j];
    for(PetscInt gi=stokesVariables->GI_V_start; gi<(stokesVariables->N_V_U_total_loc + stokesVariables->GI_V_start); gi++) 
    {
      if (gi == gi_dir){
        j++;
        if (j<stokesVariables->N_V_U_diric_loc)
          gi_dir = stokesVariables->GI_V_U_bound_dir[j]; }
      else
        stokesVariables->GI_V_U_non_dir[i] = gi, i++;
    }
  }
  else 
  {
    i = 0;
    for(PetscInt gi=stokesVariables->GI_V_start; gi<(stokesVariables->N_V_U_total_loc + stokesVariables->GI_V_start); gi++) {
      stokesVariables->GI_V_U_non_dir[i] = gi, i++; }
  }
  
  ierr = PetscMalloc1(stokesVariables->N_V_P_total_loc - stokesVariables->N_V_P_diric_loc, &stokesVariables->GI_V_P_non_dir); CHKERRQ(ierr);
  if (stokesVariables->N_V_P_diric_loc != 0)
  {
    i = 0, j = 0, gi_dir = stokesVariables->GI_V_P_bound_dir[j];
    for(PetscInt gi=stokesVariables->GI_V_P_start; gi<(stokesVariables->N_V_P_total_loc + stokesVariables->GI_V_P_start); gi++) 
    {
      if (gi == gi_dir){
        j++; 
        if (j<stokesVariables->N_V_P_diric_loc)
          gi_dir = stokesVariables->GI_V_P_bound_dir[j]; }
      else
        stokesVariables->GI_V_P_non_dir[i] = gi, i++;
    }
  }
  else 
  {
    i = 0;
    for(PetscInt gi=stokesVariables->GI_V_P_start; gi<(stokesVariables->N_V_P_total_loc + stokesVariables->GI_V_P_start); gi++) {
      stokesVariables->GI_V_P_non_dir[i] = gi, i++; }
  }

  return 0;
}


PetscErrorCode StokesVariablesPrintInfo(System_Struct *system, StokesVariables_Struct *stokesVariables, PetscBool boolStokesVariables) 
{  
  if (boolStokesVariables)
  {
    PetscPrintf(PETSC_COMM_WORLD,"  Stokes variables Info:\n");
    PetscPrintf(PETSC_COMM_WORLD,"  Number of variables %d: %d U & %d P \n", stokesVariables->Nv, stokesVariables->Nv_U, stokesVariables->Nv_P);
    PetscPrintf(PETSC_COMM_WORLD,"  Total global variables: %d\n", stokesVariables->N_V_total_glo);
    PetscPrintf(PETSC_COMM_WORLD,"  Total global U variables: %d\n", stokesVariables->N_V_U_total_glo);
    PetscPrintf(PETSC_COMM_WORLD,"  Total global P variables: %d\n", stokesVariables->N_V_P_total_glo);

    PetscPrintf(PETSC_COMM_WORLD,"  Total local variables:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d", system->rank, stokesVariables->N_V_total_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local internal variables:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d", system->rank, stokesVariables->N_V_inter_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local boundary variables:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d", system->rank, stokesVariables->N_V_bound_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local ghost variables:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d", system->rank, stokesVariables->N_V_ghost_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

    PetscPrintf(PETSC_COMM_WORLD,"  Total local U variables:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d", system->rank, stokesVariables->N_V_U_total_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local P variables:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d", system->rank, stokesVariables->N_V_P_total_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  
    for(PetscInt i=0; i<stokesVariables->Nv_U; i++){
      PetscPrintf(PETSC_COMM_WORLD,"  Global variables for U%d:\n", i);
      PetscIntView(stokesVariables->N_U_total_loc, stokesVariables->GI_V[i], stdout_viewer);}

    PetscPrintf(PETSC_COMM_WORLD,"  Global variables for P:\n");
    PetscIntView(stokesVariables->N_P_total_loc, stokesVariables->GI_V[stokesVariables->Nv_U], stdout_viewer);
  }

  return 0; 
}

PetscErrorCode StokesVariablesCreate(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements, StokesVariables_Struct *stokesVariables) 
{
  StokesVariablesBuild(system, nodes, stokesVariables);
  StokesVariablesReadBoundary(system, nodes, stokesVariables);
  StokesVariableseSearchForDirichlet(nodes, elements, stokesVariables);
  
  // Print Stokes Variables Info
  StokesVariablesPrintInfo(system, stokesVariables, print_stokesVariables);
  return 0;
}



#endif