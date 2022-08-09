#ifndef System_Functions_C
#define System_Functions_C

void SystemInitialize(System_Struct *system)
{
  /* Initialized the variables in system */
  system->FLImp_dt      = 0;
  system->FLImp_eps     = 0;

  system->FLImpExp_rTol = 1e-1;

  system->SDReload = PETSC_FALSE;   
  system->SDReload_time = 0.0;      
  system->SDReload_saved_step = 0; 

}

#endif 
