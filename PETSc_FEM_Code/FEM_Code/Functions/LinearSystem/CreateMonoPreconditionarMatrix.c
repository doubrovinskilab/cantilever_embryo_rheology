
#ifndef CreateMonolithicPreconditionarMatrix_C
#define CreateMonolithicPreconditionarMatrix_C

PetscErrorCode ComputeNonDirichletVariables(System_Struct *system, StokesVariables_Struct *stokesVariables,
  PetscInt *mom_nonDir_size_glo, PetscInt *mom_nonDir_size_loc, PetscInt *cont_nonDir_size_glo, PetscInt *cont_nonDir_size_loc)
{
  PetscErrorCode  ierr;
  PetscInt i;
  const PetscInt size = system->size;
  const PetscInt rank = system->rank;

  *mom_nonDir_size_loc =  stokesVariables->N_V_U_total_loc - stokesVariables->N_V_U_diric_loc; 
  *cont_nonDir_size_loc = stokesVariables->N_V_P_total_loc - stokesVariables->N_V_P_diric_loc; 
 
  if (size == 1) 
    *mom_nonDir_size_glo = *mom_nonDir_size_loc, *cont_nonDir_size_glo = *cont_nonDir_size_loc;
  else
  {
    *mom_nonDir_size_glo = 0, *cont_nonDir_size_glo = 0;

    
    for (i=0; i< size; i++)
      if (rank == i) { 
        *mom_nonDir_size_glo += *mom_nonDir_size_loc;
        ierr = MPI_Bcast( mom_nonDir_size_glo, 1, MPIU_INT, i, PETSC_COMM_WORLD ); CHKERRMPI(ierr);}
      else {          
        ierr = MPI_Bcast( mom_nonDir_size_glo, 1, MPIU_INT, i, PETSC_COMM_WORLD ); CHKERRMPI(ierr);}
      
    
    for (i=0; i< size; i++)
      if (rank == i) { 
        *cont_nonDir_size_glo += *cont_nonDir_size_loc;
        ierr = MPI_Bcast( cont_nonDir_size_glo, 1, MPIU_INT, i, PETSC_COMM_WORLD ); CHKERRMPI(ierr);}
      else {          
        ierr = MPI_Bcast( cont_nonDir_size_glo, 1, MPIU_INT, i, PETSC_COMM_WORLD ); CHKERRMPI(ierr); }
  }
  return 0;
}

PetscErrorCode ExtractNonDirRowsMom(System_Struct *system, StokesVariables_Struct *stokesVariables, 
  const PetscInt mom_nonDir_size_loc, const PetscInt cont_nonDir_size_loc, PetscInt *indices, PetscInt *mom_indices)
{
  PetscInt i, j, k, Li, Gi;

  j = 0, k = 0;
  if (stokesVariables->N_V_U_diric_loc != 0)
    Li = stokesVariables->GI_V_U_bound_dir[j] - stokesVariables->GI_V_start; 
  else
    Li = -1; 

  for (i=0; i < stokesVariables->N_V_U_total_loc; i++) 
  { 
    if (i==Li) 
    {
      j++; 
      if (j < stokesVariables->N_V_U_diric_loc)
        Li = stokesVariables->GI_V_U_bound_dir[j] - stokesVariables->GI_V_start;
    }
    else 
    { 
      Gi = i + stokesVariables->GI_V_start;
      mom_indices[k] = Gi;
      indices[k] = Gi;
      k++; 
    }
  }

  if (k != mom_nonDir_size_loc){
    PetscErrorPrintf("Error k = %d and should be %d on rank %d\n", k, mom_nonDir_size_loc, system->rank);
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error ExtractNonDirRowsMom"); }

  
  return 0;
}

PetscErrorCode ExtractNonDirRowsCont(System_Struct *system, StokesVariables_Struct *stokesVariables, 
  const PetscInt mom_nonDir_size_loc, const PetscInt cont_nonDir_size_loc, PetscInt *indices, PetscInt *cont_indices)
{
  PetscInt i, j, k, Li, Gi;
  
  j = 0, k = 0;
  if (stokesVariables->N_V_P_diric_loc != 0)
    Li = stokesVariables->GI_V_P_bound_dir[j] - stokesVariables->GI_V_P_start;
  else
    Li = -1;

  for (i=0; i < stokesVariables->N_V_P_total_loc; i++) 
  { 
    if (i==Li)
    {
      j++; 
      if (j < stokesVariables->N_V_P_diric_loc)
        Li = stokesVariables->GI_V_P_bound_dir[j] - stokesVariables->GI_V_P_start;
    }
    else
    {
      Gi = i + stokesVariables->GI_V_P_start;
      cont_indices[k] = Gi;
      indices[k+mom_nonDir_size_loc] = Gi;
      k++;
    } 
  }

  if (k != cont_nonDir_size_loc){
    PetscErrorPrintf("Error k = %d and should be %d on rank %d \n", k, cont_nonDir_size_loc, system->rank);
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error ExtractNonDirRowsCont"); }

  return 0;
}


PetscErrorCode ConvertStokesFullLSToReducedMonoLS(System_Struct *system, Nodes_Struct *nodes, 
  Elements_Struct *elements, StokesVariables_Struct *stokesVariables,
  FullLinearSystem_Stokes_Struct *stokesFullLS,
  ReducedMonoLinearSystem_Stokes_Struct *stokesMonoLS)
{
  PetscErrorCode  ierr;
  PetscLogDouble Clock1, Clock2;
  PetscTime(&Clock1);
  
  PetscInt mom_nonDir_size_glo, cont_nonDir_size_glo;  
  PetscInt mom_nonDir_size_loc, cont_nonDir_size_loc;

  PetscInt  *indices;       
  PetscInt  *mom_indices;   
  PetscInt  *cont_indices;  
    
  PetscPrintf(PETSC_COMM_WORLD,"Converting the Full Matrix A to a reduced monolithic matrix Ar\n");
  PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------\n"); 

  ComputeNonDirichletVariables(system, stokesVariables, &mom_nonDir_size_glo, &mom_nonDir_size_loc, &cont_nonDir_size_glo, &cont_nonDir_size_loc);
  
  PetscMalloc1(mom_nonDir_size_loc + cont_nonDir_size_loc, &indices);
  PetscMalloc2(mom_nonDir_size_loc, &mom_indices, cont_nonDir_size_loc, &cont_indices);
 
  ExtractNonDirRowsMom(system, stokesVariables, mom_nonDir_size_loc, cont_nonDir_size_loc, indices, mom_indices);
  ierr = PetscBarrier(NULL); CHKERRQ(ierr);

  ExtractNonDirRowsCont(system, stokesVariables, mom_nonDir_size_loc, cont_nonDir_size_loc, indices, cont_indices);
  ierr = PetscBarrier(NULL); CHKERRQ(ierr);
    
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, mom_nonDir_size_loc+cont_nonDir_size_loc, indices, PETSC_COPY_VALUES, &stokesMonoLS->is_nonDir_FullA);CHKERRQ(ierr);
  ierr = PetscBarrier( (PetscObject) stokesMonoLS->is_nonDir_FullA); CHKERRQ(ierr);

  ierr = ISCreateGeneral(PETSC_COMM_WORLD, mom_nonDir_size_loc, mom_indices, PETSC_COPY_VALUES, &stokesMonoLS->is_nonDir_FullA_mom);CHKERRQ(ierr);
  ierr = PetscBarrier( (PetscObject) stokesMonoLS->is_nonDir_FullA_mom); CHKERRQ(ierr);

  ierr = ISCreateGeneral(PETSC_COMM_WORLD, cont_nonDir_size_loc, cont_indices, PETSC_COPY_VALUES, &stokesMonoLS->is_nonDir_FullA_cont);CHKERRQ(ierr);
  ierr = PetscBarrier( (PetscObject) stokesMonoLS->is_nonDir_FullA_cont); CHKERRQ(ierr);
  
  PetscPrintf(PETSC_COMM_WORLD,"  Extracting the global monolithic matrix from Full matrix\n");
  
  ierr = MatCreateSubMatrix(stokesFullLS->A, stokesMonoLS->is_nonDir_FullA, stokesMonoLS->is_nonDir_FullA, 
    MAT_INITIAL_MATRIX, &stokesMonoLS->Ar); CHKERRQ(ierr);
    
  if (check_symmetry_Ar){
    PetscBool  isSymmetric;  
    PetscPrintf(PETSC_COMM_WORLD,"  Checking if symmetric.\n");
    ierr = MatIsSymmetric(stokesMonoLS->Ar, EPS, &isSymmetric); CHKERRQ(ierr);
    if (isSymmetric)  PetscPrintf(PETSC_COMM_WORLD,"    Matrix A is symmetric.\n");
    else PetscPrintf(PETSC_COMM_WORLD,"    WARNING!! Matrix A is non-symmetric.\n");
  }

  PetscPrintf(PETSC_COMM_WORLD,"  Creating the RHS vector with Dirichlet variables\n");

  PetscScalar  *zero_indices;    
  PetscMalloc1(mom_nonDir_size_loc + cont_nonDir_size_loc, &zero_indices);
  for (PetscInt i=0; i<mom_nonDir_size_loc + cont_nonDir_size_loc; i++) zero_indices[i] = 0.0;

  Vec c;
  ierr = VecDuplicate(stokesFullLS->b, &c); CHKERRQ(ierr);  
  ierr = VecCopy(stokesFullLS->b,  c); CHKERRQ(ierr);

  ierr = VecSetValues(c, mom_nonDir_size_loc+cont_nonDir_size_loc, indices, zero_indices,INSERT_VALUES); 

  ierr = VecDuplicate(stokesFullLS->b, &stokesMonoLS->b_dir); CHKERRQ(ierr);  
  ierr = MatMult(stokesFullLS->A , c, stokesMonoLS->b_dir); 

  ierr = VecScale(stokesMonoLS->b_dir, -1.0); 

  ierr = VecDestroy(&c);
  
  ierr = PetscFree(zero_indices); CHKERRQ(ierr);
  ierr = PetscFree(indices); CHKERRQ(ierr);
  ierr = PetscFree2(mom_indices, cont_indices);CHKERRQ(ierr);

  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"  Computational time: %lf seconds.\n\n", (Clock2 - Clock1) );
  
  if (print_mat_Ar){
    PetscPrintf(PETSC_COMM_WORLD,"The Monolithic reduced matrix Ar:\n"); 
    MatView(stokesMonoLS->Ar, stdout_viewer); }

  if (print_vec_b_dir){
    PetscPrintf(PETSC_COMM_WORLD,"The Monolithic reduced vector b_dir:\n"); 
    MatView(stokesMonoLS->Ar, stdout_viewer); }

  return 0;
}

PetscErrorCode CreateMonolithicPreconditionarMatrix(System_Struct *system, Nodes_Struct *nodes, 
  Elements_Struct *elements, StokesVariables_Struct *stokesVariables,
  FullLinearSystem_Stokes_Struct *stokesFullLS,
  ReducedMonoLinearSystem_Stokes_Struct *stokesMonoLS)
{
  PetscErrorCode ierr;
  PetscLogDouble Clock1, Clock2;
  PetscTime(&Clock1);

  const PetscInt Dim  = system->Dim;
  const PetscInt Nv_U = stokesVariables->Nv_U;     

  const FLUID_ELEMENT_TYPE P_ET = elements->P_elem_type; 
  const FLUID_ELEMENT_TYPE U_ET = elements->U_elem_type; 

  const PetscInt P_S_inter = elements->P_S_inter;  
  
  const PetscInt U_S_inter = elements->U_S_inter;  
  

  PetscInt e, i, n;

  PetscScalar **Muek;  
  PetscScalar **Mpek;  
  PetscScalar *d1, *d2;

  PetscInt  *GI_U_Int, *GI_P_Int;

  PetscInt mom_nonDir_size_glo, cont_nonDir_size_glo; 
  PetscInt mom_nonDir_size_loc, cont_nonDir_size_loc; 

  PetscInt first_row, last_row;

  Mat Mass, Mass_Ar; 
  Vec Diag_Ar, Diag_Ar2, B;
  Vec Diag_Mass, Diag_Mass2;
 
  MassFunction Calc_Muek = DefineMassFunction(Dim, U_ET);
  MassFunction Calc_Mpek = DefineMassFunction(Dim, P_ET);
  
  PetscMalloc1(U_S_inter, &Muek);
  for(i = 0; i < U_S_inter; ++i) PetscMalloc1(U_S_inter, &Muek[i]);

  PetscMalloc1(P_S_inter, &Mpek);
  for(i = 0; i < P_S_inter; ++i) PetscMalloc1(P_S_inter, &Mpek[i]);

  PetscMalloc2(U_S_inter, &GI_U_Int, P_S_inter, &GI_P_Int);

  ComputeNonDirichletVariables(system, stokesVariables, &mom_nonDir_size_glo, &mom_nonDir_size_loc, &cont_nonDir_size_glo, &cont_nonDir_size_loc);
  

  CreateFullPreDeterminedStokesGlobalMatrix(&Mass, system, nodes, elements, stokesVariables);

  for (e = 0; e < elements->N_inter_loc; e++)
  {
    Calc_Muek(e, nodes->x, nodes->y, nodes->z, elements->inter, Muek);
    Calc_Mpek(e, nodes->x, nodes->y, nodes->z, elements->inter, Mpek);

    
    for (n = 0; n < Nv_U; n++) 
    {
      for (i=0; i<U_S_inter;++i)  GI_U_Int[i] = stokesVariables->GI_V[n][elements->inter[e][i]]; 

      for(i=0; i< U_S_inter; i++) 
      {
        const PetscInt li   = elements->inter[e][i];   
        const PetscInt gn_u = stokesVariables->GI_V[n][li];  
        MatSetValues(Mass, 1,  &gn_u, U_S_inter, GI_U_Int, Muek[i], ADD_VALUES);  
      }
    }
    
    for (i=0; i<P_S_inter;++i)  GI_P_Int[i] = stokesVariables->GI_V[Nv_U][elements->inter[e][i]]; 

    for(i=0; i< P_S_inter; i++) 
    {
      const PetscInt li_p = nodes->LI_UtoP[elements->inter[e][i]];
      const PetscInt gn_p = stokesVariables->GI_V[n][li_p]; 
      MatSetValues(Mass, 1,  &gn_p, P_S_inter, GI_P_Int, Mpek[i], ADD_VALUES);  
    }
  }

  ierr = MatAssemblyBegin(Mass, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(Mass, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  MatZeroRows(Mass, stokesVariables->N_V_U_diric_loc, stokesVariables->GI_V_U_bound_dir, 0.0, NULL, NULL); 
  MatZeroRows(Mass, stokesVariables->N_V_P_diric_loc, stokesVariables->GI_V_P_bound_dir, 0.0, NULL, NULL); 

  ierr = MatCreateSubMatrix(Mass, stokesMonoLS->is_nonDir_FullA, stokesMonoLS->is_nonDir_FullA, MAT_INITIAL_MATRIX, &Mass_Ar); CHKERRQ(ierr);

  ierr = MatCreateAIJ(PETSC_COMM_WORLD, mom_nonDir_size_loc+cont_nonDir_size_loc, mom_nonDir_size_loc+cont_nonDir_size_loc,
          mom_nonDir_size_glo+cont_nonDir_size_glo, mom_nonDir_size_glo+cont_nonDir_size_glo, 1, NULL, 0, NULL, &stokesMonoLS->Pr); CHKERRQ(ierr);
  ierr = MatSetFromOptions(stokesMonoLS->Pr); CHKERRQ(ierr);
  ierr = MatSetUp(stokesMonoLS->Pr);  CHKERRQ(ierr);

  ierr = VecGetSubVector(stokesFullLS->b, stokesMonoLS->is_nonDir_FullA, &B);      CHKERRQ(ierr); 
  ierr = VecDuplicate(B, &Diag_Ar);                           CHKERRQ(ierr); 
  ierr = VecDuplicate(B, &Diag_Mass);                         CHKERRQ(ierr); 
  ierr = VecRestoreSubVector(stokesFullLS->b, stokesMonoLS->is_nonDir_FullA, &B);  CHKERRQ(ierr); 

  MatGetDiagonal(stokesMonoLS->Ar, Diag_Ar);
  VecDuplicate(Diag_Ar, &Diag_Ar2);
  VecCopy(Diag_Ar, Diag_Ar2);
  VecDestroy(&Diag_Ar);
  VecGetArray(Diag_Ar2, &d1);
  MatGetDiagonal(Mass_Ar, Diag_Mass);
  VecDuplicate(Diag_Mass, &Diag_Mass2);
  VecCopy(Diag_Mass, Diag_Mass2);
  VecDestroy(&Diag_Mass);
  VecGetArray(Diag_Mass2, &d2);
  MatGetOwnershipRange(stokesMonoLS->Pr, &first_row, &last_row);
  for (PetscInt r=first_row; r<last_row; r++){
    i = r - first_row;
    if (d1[i] != 0)       MatSetValue(stokesMonoLS->Pr, r, r, d1[i], INSERT_VALUES); 
    else                  MatSetValue(stokesMonoLS->Pr, r, r, d2[i], INSERT_VALUES); 
  }
  VecRestoreArray(Diag_Ar2, &d1); VecRestoreArray(Diag_Mass2, &d2);

  ierr = MatAssemblyBegin(stokesMonoLS->Pr, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);
  ierr = MatAssemblyEnd(stokesMonoLS->Pr, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr);

  VecDestroy(&Diag_Ar2);
  VecDestroy(&Diag_Mass2);

  
  for(i = 0; i < U_S_inter; ++i) PetscFree(Muek[i]);
  PetscFree(Muek);

  for(i = 0; i < P_S_inter; ++i) PetscFree(Mpek[i]);
  PetscFree(Mpek);

  PetscFree2(GI_U_Int, GI_P_Int);

  MatDestroy(&Mass_Ar);
  MatDestroy(&Mass);
  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"  Computational time: %lf seconds.\n\n", (Clock2 - Clock1) ); 

  if (print_mat_Pr){
    PetscPrintf(PETSC_COMM_WORLD,"The Monolithic reduced preconditioner matrix Pr:\n"); 
    MatView(stokesMonoLS->Pr, stdout_viewer); }

  return 0;
}
#endif