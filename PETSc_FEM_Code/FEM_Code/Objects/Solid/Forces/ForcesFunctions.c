#ifndef Solid_Forces_Functions_C
#define Solid_Forces_Functions_C

#include "BoundaryForcesFunctions.c"
#include "LinearElasticityForcesFunctions.c"
#include "SoftForcesFunctions.c"

void ForcesReset(Solid_Struct *solid){ 
  for(PetscInt j=0; j< solid->N_total_pts; j++)
    solid->forces.Fx[j] = 0.0, solid->forces.Fy[j]=0.0, solid->forces.Fz[j] = 0.0;}

PetscErrorCode ForcesFree(Solid_Struct *solid) {
  PetscPrintf(PETSC_COMM_WORLD,"    Solid Forces memory is being freed\n");
  PetscFree3(solid->forces.Fx, solid->forces.Fy, solid->forces.Fz);
  if (solid->forces.BC_FS.N>0) 
    PetscFree5(solid->forces.BC_FS.p, solid->forces.BC_FS.g, solid->forces.BC_FS.x, solid->forces.BC_FS.y, solid->forces.BC_FS.z);
  if (solid->forces.BC_FSE.N>0) 
    PetscFree3(solid->forces.BC_FSE.p, solid->forces.BC_FSE.g, solid->forces.BC_FSE.theta);
  LinearElasticityFree(solid);
  return 0; }


PetscErrorCode ForcesSumAll(Solid_Struct *solid){ 
  PetscScalar Sum_Fx = 0, Sum_Fy = 0, Sum_Fz = 0;
  for(PetscInt j=0; j< solid->N_total_pts; j++)
    Sum_Fx += solid->forces.Fx[j], Sum_Fy += solid->forces.Fy[j], Sum_Fz += solid->forces.Fz[j];
  PetscPrintf(PETSC_COMM_WORLD,"    Sum of all the forces is Total_F: %lf, %lf, %lf\n", Sum_Fx, Sum_Fy, Sum_Fz);
  if(PetscIsNanReal(Sum_Fx) || PetscIsNanReal(Sum_Fy) || PetscIsNanReal(Sum_Fz)){
    PetscErrorPrintf("Error!! Forces are Nan.\n");
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error!! Forces are Nan.");  }
  return 0; }

PetscErrorCode ForcesMallocPointArrays(Solid_Struct *solid) {
  if(solid->N_total_pts <= 0) 
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error solid->N_total_pts"); 
  
  PetscMalloc3(solid->N_total_pts, &solid->forces.Fx, solid->N_total_pts, &solid->forces.Fy, solid->N_total_pts, &solid->forces.Fz);
  return 0; }
#endif

#include "ForcesActiveFunctions.c"
