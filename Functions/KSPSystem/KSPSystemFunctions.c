#ifndef KSPSystemFunctions_H
#define KSPSystemFunctions_H

PetscErrorCode CreateKSPForStokesFullLS(FullLinearSystem_Stokes_Struct *stokesFullLS, KSP *ksp, PC *pc)
{
  PetscErrorCode ierr;
  PetscReal rtol, abstol, dtol;
  PetscInt maxits;
  PetscLogDouble Clock1, Clock2;
  KSPType ksptype;
  PCType pctype;
  PetscTime(&Clock1);
  PetscPrintf(PETSC_COMM_WORLD,"Setting up the full implicit Linear-System\n"); 
  PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------------\n"); 
  ierr = KSPCreate(PETSC_COMM_WORLD, ksp); CHKERRQ(ierr); 
  ierr = KSPSetOperators(*ksp, stokesFullLS->A, stokesFullLS->A); CHKERRQ(ierr); 

  ierr = KSPSetFromOptions(*ksp); CHKERRQ(ierr);
  ierr = KSPGetPC(*ksp, pc); CHKERRQ(ierr);
  
  ierr = KSPSetUp(*ksp); CHKERRQ(ierr);

  ierr = KSPGetTolerances(*ksp, &rtol, &abstol, &dtol, &maxits); CHKERRQ(ierr);
  ierr = KSPGetType(*ksp,  &ksptype); CHKERRQ(ierr);
  ierr = PCGetType(*pc, &pctype); CHKERRQ(ierr);
  
  PetscPrintf(PETSC_COMM_WORLD,"  KSP type = %s\n", ksptype);
  PetscPrintf(PETSC_COMM_WORLD,"  PC type = %s\n", pctype);
  PetscPrintf(PETSC_COMM_WORLD,"  Max iterations = %3D\n", maxits);
  PetscPrintf(PETSC_COMM_WORLD,"  Relative tolerance iterations = %g\n", (double) rtol);
  PetscPrintf(PETSC_COMM_WORLD,"  Absolute tolerance iterations = %g\n", (double) abstol);
  PetscPrintf(PETSC_COMM_WORLD,"  Divergence tolerance iterations = %g\n", (double) dtol);
  ierr = PetscBarrier(NULL); CHKERRQ(ierr);

  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"  Computational time: %lf seconds\n\n", (Clock2 - Clock1) );

  return 0;
}

PetscErrorCode CreateKSPForStokesReducedMonoLS(ReducedMonoLinearSystem_Stokes_Struct *stokesMonoLS, KSP *ksp, PC *pc, PetscBool boolMonoPre)
{
  PetscErrorCode ierr;
  PetscReal rtol, abstol, dtol;
  PetscInt maxits;
  PetscLogDouble Clock1, Clock2;
  KSPType ksptype;
  PCType pctype;
  PetscTime(&Clock1);

  
  PetscPrintf(PETSC_COMM_WORLD,"Setting up the reduced mono Linear-System...\n"); 
  PetscPrintf(PETSC_COMM_WORLD,"------------------------------------------------\n"); 
  ierr = KSPCreate(PETSC_COMM_WORLD,ksp); CHKERRQ(ierr);
  if (boolMonoPre){ ierr = KSPSetOperators(*ksp, stokesMonoLS->Ar, stokesMonoLS->Pr); CHKERRQ(ierr);  }
  else            { ierr = KSPSetOperators(*ksp, stokesMonoLS->Ar, stokesMonoLS->Ar); CHKERRQ(ierr);  }
  ierr = KSPSetFromOptions(*ksp); CHKERRQ(ierr);
  ierr = KSPGetPC(*ksp, pc); CHKERRQ(ierr);
  ierr = KSPSetUp(*ksp);

  ierr = KSPGetTolerances(*ksp, &rtol, &abstol, &dtol, &maxits); CHKERRQ(ierr);
  ierr = KSPGetType(*ksp, &ksptype); CHKERRQ(ierr);
  ierr = PCGetType(*pc,&pctype); CHKERRQ(ierr);
  
  PetscPrintf(PETSC_COMM_WORLD,"  KSP type = %s\n", ksptype);
  PetscPrintf(PETSC_COMM_WORLD,"  PC type = %s\n", pctype);
  PetscPrintf(PETSC_COMM_WORLD,"  Max iterations = %3D\n", maxits);
  PetscPrintf(PETSC_COMM_WORLD,"  Relative tolerance iterations = %g\n", (double) rtol);
  PetscPrintf(PETSC_COMM_WORLD,"  Absolute tolerance iterations = %g\n", (double) abstol);
  PetscPrintf(PETSC_COMM_WORLD,"  Divergence tolerance iterations = %g\n", (double) dtol);
  ierr = PetscBarrier((PetscObject) stokesMonoLS->Ar);CHKERRQ(ierr);

  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"  Computational time: %lf seconds\n\n", (Clock2 - Clock1) ); 

  return 0;
}


PetscErrorCode CreateKSPForStokesLS(System_Struct *system, FullLinearSystem_Stokes_Struct *stokesFullLS, 
  ReducedMonoLinearSystem_Stokes_Struct *stokesMonoLS, KSP *ksp, PC *pc)
{
  if(system->FLSolver_Method == DIRECT)
    CreateKSPForStokesFullLS(stokesFullLS, ksp, pc);
  else if(system->FLSolver_Method == MONO)
    CreateKSPForStokesReducedMonoLS(stokesMonoLS, ksp, pc, PETSC_FALSE);
  else if(system->FLSolver_Method == MONO_PRE)
    CreateKSPForStokesReducedMonoLS(stokesMonoLS, ksp, pc, PETSC_TRUE);
  return 0;
}

PetscErrorCode SolveKSPForStokesFullLS(FullLinearSystem_Stokes_Struct *stokesFullLS, KSP *ksp, PC *pc)
{ 
  PetscErrorCode ierr;
  PetscReal norm_L2; 

  PetscLogDouble Clock1, Clock2;
  PetscInt i;
  Vec error; 

  PetscTime(&Clock1);

  ierr = VecDuplicate(stokesFullLS->b, &error);  CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"Solving the Full Stokes Linear-System implicitly...\n"); 

  KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE);
  KSPSetReusePreconditioner(*ksp, PETSC_TRUE);
  
  ierr = KSPSetOperators(*ksp, stokesFullLS->A, stokesFullLS->A); CHKERRQ(ierr);
  ierr = KSPSolve(*ksp, stokesFullLS->b, stokesFullLS->x); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"  System solved\n"); 
  
  MatMult(stokesFullLS->A, stokesFullLS->x, error); 
  VecAXPY(error, -1.0, stokesFullLS->b);            
  VecNorm(error, NORM_2, &norm_L2);           
  KSPGetIterationNumber(*ksp, &i);

  ierr = VecGhostUpdateBegin(stokesFullLS->x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(stokesFullLS->x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = VecDestroy(&error); CHKERRQ(ierr);
  PetscTime(&Clock2);

  PetscPrintf(PETSC_COMM_WORLD,"    Number of iterations = %3D\n",i);
  
  PetscPrintf(PETSC_COMM_WORLD,"    Error-L2 = %g\n", (double) norm_L2);
  PetscPrintf(PETSC_COMM_WORLD,"    Computational time: %lf seconds\n", (Clock2 - Clock1) );

  return 0;
}

PetscErrorCode SolveKSPForStokesReducedMonoLS(FullLinearSystem_Stokes_Struct *stokesFullLS, 
  ReducedMonoLinearSystem_Stokes_Struct *stokesMonoLS, KSP *ksp, PC *pc, PetscBool boolMonoPre)
{
  PetscErrorCode ierr;
  PetscReal norm_L2; 

  PetscLogDouble Clock1, Clock2;
  PetscInt i;

  Vec error;    
  Vec X;        
  Vec B;        
  Vec b_update; 
  PetscTime(&Clock1);

  ierr = VecDuplicate(stokesFullLS->b, &b_update); CHKERRQ(ierr);       
  ierr = VecWAXPY(b_update, 1.0, stokesFullLS->b, stokesMonoLS->b_dir); CHKERRQ(ierr); 

  ierr = VecGetSubVector(stokesFullLS->x, stokesMonoLS->is_nonDir_FullA, &X); CHKERRQ(ierr); 
  ierr = VecGetSubVector(b_update,        stokesMonoLS->is_nonDir_FullA, &B); CHKERRQ(ierr); 
  ierr = VecDuplicate(X, &error);                                             CHKERRQ(ierr); 

  if(print_vec_br) {
    VecGetSize(B,&i);
    PetscPrintf(PETSC_COMM_WORLD,"The Monolithic reduced vector br of size %d:\n", i); 
    VecView(B, stdout_viewer); }
  
  PetscPrintf(PETSC_COMM_WORLD,"  Solving the Stokes Reduced Mono Linear-System implicitly...\n"); 

  KSPSetInitialGuessNonzero(*ksp, PETSC_TRUE);
  KSPSetReusePreconditioner(*ksp, PETSC_TRUE);
  if (boolMonoPre){ ierr = KSPSetOperators(*ksp, stokesMonoLS->Ar, stokesMonoLS->Pr); CHKERRQ(ierr);  }
  else            { ierr = KSPSetOperators(*ksp, stokesMonoLS->Ar, stokesMonoLS->Ar); CHKERRQ(ierr);  }

  ierr = KSPSolve(*ksp, B, X); CHKERRQ(ierr);

  PetscPrintf(PETSC_COMM_WORLD,"  System solved\n"); ; 

  MatMult(stokesMonoLS->Ar, X, error);         
  VecAXPY(error,  -1.0,      B);               
  
  VecNorm(error,  NORM_2,        &norm_L2);    
  
  KSPGetIterationNumber(*ksp, &i);
  
  ierr = VecRestoreSubVector(stokesFullLS->x, stokesMonoLS->is_nonDir_FullA, &X); CHKERRQ(ierr);
  ierr = VecRestoreSubVector(b_update, stokesMonoLS->is_nonDir_FullA, &B); CHKERRQ(ierr);
  
  ierr = VecGhostUpdateBegin(stokesFullLS->x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);
  ierr = VecGhostUpdateEnd(stokesFullLS->x,INSERT_VALUES,SCATTER_FORWARD); CHKERRQ(ierr);

  ierr = VecDestroy(&b_update);  CHKERRQ(ierr);
  ierr = VecDestroy(&error);  CHKERRQ(ierr);
  ierr = VecDestroy(&X);  CHKERRQ(ierr);
  ierr = VecDestroy(&B);  CHKERRQ(ierr);
  PetscTime(&Clock2);

  PetscPrintf(PETSC_COMM_WORLD,"    Number of iterations = %3D\n",i);
  
  PetscPrintf(PETSC_COMM_WORLD,"    Error-L2 = %g\n", (double) norm_L2);
  PetscPrintf(PETSC_COMM_WORLD,"    Computational time: %lf seconds\n", (Clock2 - Clock1) );

  return 0;
}

PetscErrorCode SolveKSPForStokesLS(System_Struct *system, FullLinearSystem_Stokes_Struct *stokesFullLS, 
  ReducedMonoLinearSystem_Stokes_Struct *stokesMonoLS, KSP *ksp, PC *pc)
{
  if(system->FLSolver_Method == DIRECT)
    SolveKSPForStokesFullLS(stokesFullLS, ksp, pc);
  else if(system->FLSolver_Method == MONO)
    SolveKSPForStokesReducedMonoLS(stokesFullLS, stokesMonoLS, ksp, pc, PETSC_FALSE);
  else if(system->FLSolver_Method == MONO_PRE)
    SolveKSPForStokesReducedMonoLS(stokesFullLS, stokesMonoLS, ksp, pc, PETSC_TRUE);
  return 0;
}


#endif