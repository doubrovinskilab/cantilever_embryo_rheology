#ifndef AssembleStokesFullLinearSystem_C
#define AssembleStokesFullLinearSystem_C

PetscErrorCode AssembleStokesFullLS(System_Struct *system, Nodes_Struct *nodes, 
  Elements_Struct *elements, StokesVariables_Struct *stokesVariables,
  FullLinearSystem_Stokes_Struct *stokesFullLS)
{
  PetscErrorCode ierr;
  PetscLogDouble Clock1, Clock2;
  PetscTime(&Clock1);

  const PetscMPIInt size  = system->size;
  const PetscMPIInt rank  = system->rank;

  const PetscInt Dim      = system->Dim;
  const PetscInt Nv       = stokesVariables->Nv;        
  const PetscInt Nv_U     = stokesVariables->Nv_U;      
  

  const FLUID_ELEMENT_TYPE P_ET = elements->P_elem_type; 
  const FLUID_ELEMENT_TYPE U_ET = elements->U_elem_type; 

  const PetscInt P_S_inter = elements->P_S_inter;  
  const PetscInt P_S_bound = elements->P_S_bound;  
  const PetscInt U_S_inter = elements->U_S_inter;  
  const PetscInt U_S_bound = elements->U_S_bound;  

  const PetscInt N_bound_loc = elements->N_bound_loc; 

  const PetscScalar mu  = stokesVariables->mu;        

  const PetscScalar eps = system->FLImp_eps; 

  const PetscInt N_V_total_glo = stokesVariables->N_V_total_glo; 
  const PetscInt N_V_total_loc = stokesVariables->N_V_total_loc; 
  const PetscInt N_V_ghost_loc = stokesVariables->N_V_ghost_loc; 


  PetscBool boolDirPressBoundary = PETSC_FALSE; 

  PetscInt i, j;

  PetscScalar Mu;                            

  PetscScalar **Dek;                         
  PetscScalar **Gek_x,  **Gek_y,  **Gek_z;   
  PetscScalar **Gek_xT, **Gek_yT, **Gek_zT;  
  PetscScalar  *Fek_x,   *Fek_y,   *Fek_z;   
  PetscScalar **Mpek;                        
  
  PetscScalar **Plp;                         
  PetscScalar **Urlp;                        
  PetscScalar  *Unlp;                        
  PetscScalar **Pnlp;                        
  PetscScalar  *Pnlp_r;                      

  PetscInt **GI_V_elem;   
                          

  PetscInt  *GI_P;        
  PetscInt  *GI_U;        

  PetscScalar  *VecN, *VecT1, *VecT2;  
                          
  PetscPrintf(PETSC_COMM_WORLD,"Creating the Stokes Global Matrix & Vectors with eps %lf\n",eps); 
  PetscPrintf(PETSC_COMM_WORLD,"--------------------------------------------------------------\n"); 

  PetscMalloc1(U_S_inter, &Dek);
  for(i = 0; i < U_S_inter; ++i) PetscMalloc1(U_S_inter, &Dek[i]);

  PetscMalloc3(U_S_inter, &Gek_x, U_S_inter, &Gek_y, U_S_inter, &Gek_z);
  for(i = 0; i < U_S_inter; ++i) PetscMalloc3(P_S_inter, &Gek_x[i], P_S_inter, &Gek_y[i], P_S_inter, &Gek_z[i]);

  PetscMalloc3(P_S_inter, &Gek_xT, P_S_inter, &Gek_yT, P_S_inter, &Gek_zT);
  for(i = 0; i < P_S_inter; ++i) PetscMalloc3(U_S_inter, &Gek_xT[i], U_S_inter, &Gek_yT[i], U_S_inter, &Gek_zT[i]);

  PetscMalloc3(U_S_inter, &Fek_x, U_S_inter, &Fek_y, U_S_inter, &Fek_z);

  PetscMalloc1(P_S_inter, &Mpek);
  for(i = 0; i < P_S_inter; ++i) PetscMalloc1(P_S_inter, &Mpek[i]);

  PetscMalloc1(U_S_bound, &Unlp);

  PetscMalloc1(U_S_bound, &Plp);
  for(i = 0; i < U_S_bound; ++i) PetscMalloc1(P_S_bound, &Plp[i]);

  PetscMalloc1(U_S_bound, &Urlp);
  for(i = 0; i < U_S_bound; ++i) PetscMalloc1(U_S_bound, &Urlp[i]);

  PetscMalloc1(Nv, &GI_V_elem); 
  for(i = 0; i < Nv_U; ++i)  PetscMalloc1(U_S_inter, &GI_V_elem[i]); 

  PetscMalloc1(P_S_inter, &GI_V_elem[Nv-1]);  

  PetscMalloc1(P_S_bound, &Pnlp);
  for(i = 0; i < P_S_bound; ++i) PetscMalloc1(P_S_bound, &Pnlp[i]);

  PetscMalloc1(P_S_bound, &Pnlp_r);

  PetscMalloc2(P_S_bound, &GI_P, U_S_bound, &GI_U); 

  PetscMalloc3(3, &VecN, 3, &VecT1, 3, &VecT2); 
  
  PetscPrintf(PETSC_COMM_WORLD,"  Creating the global matrices and vectors\n"); 

  CreateGlobalFullVectorWithGhostCells(    &stokesFullLS->x, N_V_total_glo, N_V_total_loc, 
    N_V_ghost_loc, stokesVariables->GI_V_ghost, size);

  CreateGlobalFullVectorWithGhostCells(    &stokesFullLS->b, N_V_total_glo, N_V_total_loc, 
    N_V_ghost_loc, stokesVariables->GI_V_ghost, size);

  CreateFullPreDeterminedStokesGlobalMatrix( &stokesFullLS->A, system, nodes, elements, 
    stokesVariables);

  DiffusionFunction Calc_Dek  = DefineDiffusionFunction(Dim, U_ET);         
  GradientFunction  Calc_Gek  = DefineGradientFunction(Dim,  U_ET, P_ET);
  GradientFunction  Calc_GekT = DefineGradientFunction(Dim,  P_ET, U_ET);
  MassFunction      Calc_Mpek = DefineMassFunction(Dim, P_ET); 

  ScaleFunction     Scale_Dek  = DefineScaleFunction(Dim, U_ET, U_ET); 
  ScaleFunction     Scale_Gek  = DefineScaleFunction(Dim, U_ET, P_ET); 
  ScaleFunction     Scale_GekT = DefineScaleFunction(Dim, P_ET, U_ET); 
  ScaleFunction     Scale_Mpek = DefineScaleFunction(Dim, P_ET, P_ET); 

  AddEqsMom  Add_Mom_Func     = DefineAddMomentumFunction(Dim);
  AddEqsCont Add_Cont_Func    = DefineAddContinuityFunction(Dim, 2);     

  CalcVisc  Calc_Visc_Func = DefineCalcViscosityFunction(stokesVariables->ViscMethod);

  Mu = mu; 

  PetscPrintf(PETSC_COMM_WORLD,"  Looping over internal elements\n");

  for (PetscInt e = 0; e < elements->N_inter_loc; e++)
  {
    
    Calc_Dek(e, nodes->x, nodes->y, nodes->z, elements->inter, Dek); 
    Calc_Gek(e, nodes->x, nodes->y, nodes->z, elements->inter, Gek_x,  Gek_y,  Gek_z); 
    Calc_GekT(e, nodes->x, nodes->y, nodes->z, elements->inter, Gek_xT, Gek_yT, Gek_zT); 
    Calc_Mpek(e, nodes->x, nodes->y, nodes->z, elements->inter, Mpek);
    
    Calc_Visc_Func(&Mu, &mu, stokesVariables, &P_S_inter, elements->inter[e], nodes->x, nodes->y, nodes->z);

    Scale_Dek(Dek,     Mu);   

    Scale_Gek(Gek_x, -1.0);   
    Scale_Gek(Gek_y, -1.0);   
    Scale_Gek(Gek_z, -1.0);   

    Scale_GekT(Gek_xT, -1.0); 
    Scale_GekT(Gek_yT, -1.0); 
    Scale_GekT(Gek_zT, -1.0); 

    Scale_Mpek(Mpek, -1*eps); 

    for (i=0; i<U_S_inter; ++i)   Fek_x[i] = 0.0, Fek_y[i] = 0.0, Fek_z[i] = 0.0;

    for (i=0; i<U_S_inter;++i)
      for (j=0; j<Nv_U;++j)
        GI_V_elem[j][i] = stokesVariables->GI_V[j][elements->inter[e][i]]; 

    j = Nv_U;
    for (i=0; i<P_S_inter;++i)
      GI_V_elem[j][i]  = stokesVariables->GI_V[j][nodes->LI_UtoP[elements->inter[e][i]]]; 

    for (i = 0; i< U_S_inter; ++i)
      Add_Mom_Func(&stokesFullLS->A, &stokesFullLS->b, U_S_inter, P_S_inter, i, GI_V_elem, Dek, Gek_x, Gek_y, Gek_z, Fek_x, Fek_y, Fek_z); 

    for (i = 0; i< P_S_inter; ++i)
      Add_Cont_Func(&stokesFullLS->A, U_S_inter, P_S_inter, i, GI_V_elem, Mpek, Gek_xT, Gek_yT, Gek_zT); 
  }

  PetscPrintf(PETSC_COMM_WORLD,"    Finished looping over internal elements\n");
  ierr = MatAssemblyBegin(stokesFullLS->A, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr); 
  ierr = MatAssemblyEnd(stokesFullLS->A, MAT_FLUSH_ASSEMBLY); CHKERRQ(ierr); 
  
  PetscPrintf(PETSC_COMM_WORLD,"  Looping over velocity boundary elements\n"); 
  for (PetscInt be=0; be < N_bound_loc; ++be)
  {
    const PetscInt g = elements->Gr_bound[be];

    for (PetscInt n = 0; n < Nv_U; n++) 
    {
      if (stokesVariables->Gr_bound_type[n][g] != DIRICHLET) 
      {

        for (i=0; i< P_S_bound; i++)
        { const PetscInt li   = elements->bound[be][i];  
          const PetscInt li_P = nodes->LI_UtoP[li];      
          GI_P[i] = stokesVariables->GI_V[Nv_U][li_P];  }
        
        for (i=0; i< U_S_bound; i++)
        { const PetscInt li   = elements->bound[be][i];  
          GI_U[i] = stokesVariables->GI_V[n][li];}

        Calc_Visc_Func(&Mu, &mu, stokesVariables, &P_S_bound, elements->bound[be], nodes->x, nodes->y, nodes->z);

        if ( stokesVariables->Gr_bound_type[n][g] == NEUMANN || stokesVariables->Gr_bound_type[n][g] == ROBIN || stokesVariables->Gr_bound_type[n][g] == ZERO_NORMAL_VELOCITY || stokesVariables->Gr_bound_type[n][g] == SLIP)
        {
          SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT , "Error only 'Dirichlet' boundary is defined for Velocity boundary for now");
        }
      }
    }
  }

  PetscPrintf(PETSC_COMM_WORLD,"    Finished looping over velocity boundary elements\n"); 
  ierr = VecAssemblyBegin(stokesFullLS->b); CHKERRQ(ierr); 
  ierr = VecAssemblyEnd(stokesFullLS->b); CHKERRQ(ierr); 
  
  PetscPrintf(PETSC_COMM_WORLD,"    Final assembly of the matrix A\n"); 
  ierr = MatAssemblyBegin(stokesFullLS->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
  ierr = MatAssemblyEnd(stokesFullLS->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
    
  PetscPrintf(PETSC_COMM_WORLD,"  Setting Dirichlet/Explicit boundary nodes in A to 0 except diag to 1\n"); 
  ierr = MatZeroRows(stokesFullLS->A, stokesVariables->N_V_U_diric_loc, stokesVariables->GI_V_U_bound_dir, 1.0, NULL, NULL); CHKERRQ(ierr); 
  ierr = MatZeroRows(stokesFullLS->A, stokesVariables->N_V_P_diric_loc, stokesVariables->GI_V_P_bound_dir, 1.0, NULL, NULL); CHKERRQ(ierr); 

  PetscPrintf(PETSC_COMM_WORLD,"  Looping over Dirichlet boundary nodes\n"); 

  for (PetscInt be=0; be < elements->N_bound_loc; ++be)
  {
    const PetscInt g  = elements->Gr_bound[be]; 

    for (PetscInt n = 0; n < Nv_U; n++)
      if (stokesVariables->Gr_bound_type[n][g] == DIRICHLET) 
        for (i = 0; i < elements->U_S_bound; i++){
          const PetscInt li = elements->bound[be][i];
          const PetscInt gn_u = stokesVariables->GI_V[n][li]; 
          ierr = VecSetValue(stokesFullLS->b, gn_u, stokesVariables->Gr_bound_value[n][g], INSERT_VALUES); CHKERRQ(ierr); 
          ierr = VecSetValue(stokesFullLS->x, gn_u, stokesVariables->Gr_bound_value[n][g], INSERT_VALUES); CHKERRQ(ierr); 
        }

    PetscInt n = Nv_U;
    if (stokesVariables->Gr_bound_type[n][g] == DIRICHLET) 
      for (i = 0; i < elements->P_S_bound; i++){
        const PetscInt li_p = nodes->LI_UtoP[elements->bound[be][i]];
        const PetscInt gn_p = stokesVariables->GI_V[n][li_p]; 
        ierr = VecSetValue(stokesFullLS->b, gn_p, stokesVariables->Gr_bound_value[n][g], INSERT_VALUES); CHKERRQ(ierr); 
        ierr = VecSetValue(stokesFullLS->x, gn_p, stokesVariables->Gr_bound_value[n][g], INSERT_VALUES); CHKERRQ(ierr); 
        }
  }

  PetscPrintf(PETSC_COMM_WORLD,"    Final assembly of the vectors x & b\n"); 
  ierr = VecAssemblyBegin(stokesFullLS->b); CHKERRQ(ierr); 
  ierr = VecAssemblyEnd(stokesFullLS->b); CHKERRQ(ierr); 
  ierr = VecAssemblyBegin(stokesFullLS->x); CHKERRQ(ierr); 
  ierr = VecAssemblyEnd(stokesFullLS->x); CHKERRQ(ierr); 

  ierr = ISCreateStride(PETSC_COMM_WORLD, stokesVariables->N_V_U_total_loc, stokesVariables->GI_V_start, 1, &stokesFullLS->is_mom); CHKERRQ(ierr);
  ierr = ISCreateStride(PETSC_COMM_WORLD, stokesVariables->N_V_P_total_loc, stokesVariables->GI_V_P_start, 1, &stokesFullLS->is_cont); CHKERRQ(ierr);
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, stokesVariables->N_V_U_total_loc - stokesVariables->N_V_U_diric_loc, 
    stokesVariables->GI_V_U_non_dir, PETSC_COPY_VALUES, &stokesFullLS->is_mom_non_dir); CHKERRQ(ierr);  
  ierr = ISCreateGeneral(PETSC_COMM_WORLD, stokesVariables->N_V_P_total_loc - stokesVariables->N_V_P_diric_loc, 
    stokesVariables->GI_V_P_non_dir, PETSC_COPY_VALUES, &stokesFullLS->is_cont_non_dir); CHKERRQ(ierr);

  boolDirPressBoundary = PETSC_FALSE;   
  for(i = 0; i < nodes->N_Gr; ++i)
    if (stokesVariables->Gr_bound_type[Nv_U][i] == DIRICHLET){
      boolDirPressBoundary = PETSC_TRUE;
      break;
    }

  if (PetscNot(boolDirPressBoundary))
  {
    PetscPrintf(PETSC_COMM_WORLD,"    Since pressure has no fixed boundary then one point becomes Dirichlet\n"); 
    PetscScalar DN[3];  
    PetscScalar DN_P;   
    PetscInt Li, gn_p;
    PetscInt rankFound;
    SearchInputAsciiFileFor3Scalar(system->SYSInputFile, "DirichletNode_Coor",    DN);
    SearchInputAsciiFileForScalar(system->SYSInputFile,  "DirichletNode_PValue", &DN_P);
    NodesSearch(system, nodes, DN, &Li,  &rankFound);
    if (rankFound == rank)
    {
      const PetscInt li_p = nodes->LI_UtoP[Li];
      gn_p = stokesVariables->GI_V[Nv_U][li_p]; 
      ierr = VecSetValue( stokesFullLS->b, gn_p, DN_P, INSERT_VALUES); CHKERRQ(ierr);
      ierr = MPI_Bcast( &gn_p, 1, MPIU_INT, rankFound, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
      printf("    The node found is: %lf, %lf, %lf\n", nodes->x[Li], nodes->y[Li], nodes->z[Li]);
    }
    else
    {
      ierr = MPI_Bcast( &gn_p, 1, MPIU_INT, rankFound, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
    }

    ierr = MatZeroRowsColumns(stokesFullLS->A, 1, &gn_p, 1.0, NULL, NULL); CHKERRQ(ierr);

    ierr = VecAssemblyBegin(stokesFullLS->b); CHKERRQ(ierr); 
    ierr = VecAssemblyEnd(stokesFullLS->b); CHKERRQ(ierr); 

    ierr = MatAssemblyBegin(stokesFullLS->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
    ierr = MatAssemblyEnd(stokesFullLS->A, MAT_FINAL_ASSEMBLY); CHKERRQ(ierr); 
    
    PetscPrintf(PETSC_COMM_WORLD,"    The pressure variable point %d on rank %d has been fixed to %lf\n", gn_p, rankFound, DN_P); 
  }

  ierr = PetscPrintf(PETSC_COMM_WORLD,"  Free memory\n"); 

  for(i = 0; i < U_S_inter; ++i) ierr = PetscFree(Dek[i]);
  ierr = PetscFree(Dek);

  for(i = 0; i < U_S_inter; ++i) ierr = PetscFree3(Gek_x[i], Gek_y[i], Gek_z[i]);
  ierr = PetscFree3(Gek_x, Gek_y, Gek_z);

  for(i = 0; i < P_S_inter; ++i) ierr = PetscFree3(Gek_xT[i], Gek_yT[i], Gek_zT[i]);
  ierr = PetscFree3(Gek_xT, Gek_yT, Gek_zT);

  for(i = 0; i < P_S_inter; ++i) ierr = PetscFree(Mpek[i]);
  ierr = PetscFree(Mpek);

  ierr = PetscFree(Unlp);

  for(i = 0; i < U_S_bound; ++i) ierr = PetscFree(Plp[i]);
  ierr = PetscFree(Plp);

  for(i = 0; i < U_S_bound; ++i) ierr = PetscFree(Urlp[i]);
  ierr = PetscFree(Urlp);

  for(i = 0; i < P_S_bound; ++i) ierr = PetscFree(Pnlp[i]);
  ierr = PetscFree(Pnlp);

  ierr = PetscFree(Pnlp_r);

  ierr = PetscFree3(Fek_x, Fek_y, Fek_z);

  for(i = 0; i < Nv; ++i) ierr = PetscFree(GI_V_elem[i]);
  ierr = PetscFree(GI_V_elem);

  ierr = PetscFree2(GI_P, GI_U);

  ierr = PetscFree3(VecN, VecT1, VecT2);

  if(print_mat_A) {
    PetscPrintf(PETSC_COMM_WORLD,"The Stokes Full Implicit A matrix is:\n"); 
    MatView(stokesFullLS->A, stdout_viewer); }

  if(print_vec_b) {
    PetscPrintf(PETSC_COMM_WORLD,"The Stokes Full Implicit b vector is:\n"); 
    VecView(stokesFullLS->b, stdout_viewer); }

  PetscTime(&Clock2);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"  Computational time: %lf seconds\n\n", (Clock2 - Clock1) ); CHKERRQ(ierr);

  return 0;
}

#endif