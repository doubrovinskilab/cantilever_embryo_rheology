#ifndef IB_Functions_H
#define IB_Functions_H

#include "SearchFunctions/SearchForBoxesWIthSolidPoints.c"
#include "TransferFunctions/TransferForcesToFluid.c"
#include "TransferFunctions/TransferVelocityToSolid.c"


PetscErrorCode IBTimeStepping(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements,
  StokesVariables_Struct *stokesVariables, FullLinearSystem_Stokes_Struct *stokesFullLS,
  ReducedMonoLinearSystem_Stokes_Struct *stokesMonoLS, Boxes_Struct *boxes, Solid_Struct *solid, 
  ForcesActive_Struct *active_forces, KSP *ksp, PC *pc)
{
  PetscErrorCode  ierr;
  system->time = 0;
  PetscScalar     dt;
  PetscScalar     dtime_save  = system->dtime_save;
  PetscScalar     end_time    = system->end_time;
  PetscScalar     time_save   = 0;
  PetscInt        time_step   = 0;
  PetscInt        saved_step  = 0;

  PetscScalar     fluid_dtime_save = 0;
  PetscScalar     fluid_time_save  = 0;
  if (write_fluid_info) fluid_dtime_save = system->dtime_save*write_fluid_N_dtime;

  PetscLogDouble Clock1, Clock2;

  ierr = VecDuplicate(stokesFullLS->b, &stokesFullLS->b0); CHKERRQ(ierr);  
  ierr = VecCopy(stokesFullLS->b, stokesFullLS->b0); CHKERRQ(ierr);

  if (system->FLSolver_Type == IMPLICIT) 
    dt = system->FLImp_dt;

  PetscPrintf(PETSC_COMM_WORLD,"\nStarting the IBM Solver:\n"); 
  PetscPrintf(PETSC_COMM_WORLD,"=====================================\n"); 
  PetscPrintf(PETSC_COMM_WORLD,"Time step size: %f sec\n", dt); 
  PetscPrintf(PETSC_COMM_WORLD,"Final time: %f\n", end_time); 
  PetscPrintf(PETSC_COMM_WORLD,"Incrimental time to save solid data: %f\n", dtime_save); 
  if PetscNot(write_fluid_info) 
    PetscPrintf(PETSC_COMM_WORLD,"Fluid data will not be saved\n"); 
  else
    PetscPrintf(PETSC_COMM_WORLD,"Incrimental time to save fluid data: %f\n", fluid_dtime_save); 
  
  system->dt = dt;

  ForcesInitialize(system, solid, active_forces);

  BaseFunction                    Base_Func  = DefineBaseFunction(system->Dim, elements->U_elem_type);

  if (system->SDReload) {
    ReadVTKSolidFileAndUpdatePoints( system,  solid);
    system->time= system->SDReload_time;  
    time_step   = (PetscInt) system->time/dt; 
    saved_step  = system->SDReload_saved_step+1; 
    time_save   = system->time+dtime_save;
    fluid_time_save = system->time+fluid_dtime_save;
    BoundaryForce_SpeedForceReload(solid, system->time);
    BoundaryForce_SpeedForceEllipseReload(solid, system->time);}

  PetscPrintf(PETSC_COMM_WORLD,"\n"); 
  for(;system->time<=end_time; system->time+=dt)
  {
    PetscTime(&Clock1);
    PetscPrintf(PETSC_COMM_WORLD,"Time step: %i, Simulation time: %f seconds\n", time_step, system->time);

    SearchForBoxWithPoints(system, nodes, elements, boxes, solid); 

    ForcesCalculate(system, solid, active_forces); 

    ForcesSumAll(solid);

    ierr = VecCopy(stokesFullLS->b0, stokesFullLS->b); CHKERRQ(ierr); 

    TransferForceToFluid(system, nodes, elements, stokesVariables, stokesFullLS, boxes, solid, &Base_Func); 

    SolveKSPForStokesLS(system, stokesFullLS, stokesMonoLS, ksp, pc);

    if (system->time >= time_save)
    { 
      VTKWriteSolidData(system, solid, saved_step);

      if (write_fluid_info && system->time >= fluid_time_save){
        VTKWriteFluidStokesVariables(system, nodes, elements, stokesVariables, &stokesFullLS->x, saved_step);
        fluid_time_save+=fluid_dtime_save; }

      time_save+= dtime_save; 
      saved_step++;
    }

    TransferVelocityFromFluidToSolidPoints(system, nodes, elements, stokesVariables, stokesFullLS, boxes, solid, &Base_Func, dt);
    
    PetscBarrier(NULL); 
    PetscTime(&Clock2);
    PetscPrintf(PETSC_COMM_WORLD,"  Time step computational time: %lf seconds\n\n", (Clock2 - Clock1) ); 
    time_step++;
  }


  return 0;
}

#endif