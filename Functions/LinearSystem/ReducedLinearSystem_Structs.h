
typedef struct {
  Mat Ar;       
  Mat Pr;       
  Vec b_dir;
  IS  is_nonDir_FullA;       
  IS  is_nonDir_FullA_mom;   
  IS  is_nonDir_FullA_cont;  
} ReducedMonoLinearSystem_Stokes_Struct;

PetscErrorCode ReducedMonoLinearSystem_Stokes_Free(System_Struct *system, ReducedMonoLinearSystem_Stokes_Struct *stokesMonoLS) 
{
  
  if(system->FLSolver_Method == MONO || system->FLSolver_Method == MONO_PRE)
  {
    MatDestroy(&stokesMonoLS->Ar);
    if (system->FLSolver_Method == MONO_PRE)
      MatDestroy(&stokesMonoLS->Pr);

    VecDestroy(&stokesMonoLS->b_dir);
    
    ISDestroy(&stokesMonoLS->is_nonDir_FullA);
    ISDestroy(&stokesMonoLS->is_nonDir_FullA_mom);
    ISDestroy(&stokesMonoLS->is_nonDir_FullA_cont);
  }

  return 0; 
}