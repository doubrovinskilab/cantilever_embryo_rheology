#ifndef Solid_Forces_LinearElasticity_Functions_C
#define Solid_Forces_LinearElasticity_Functions_C

PetscErrorCode LinearElasticity(System_Struct *system, Solid_Struct *solid)
{ 
  PetscScalar nVx=0, nVy=0, nVz=0, L, Ef;
  const PetscScalar K   = solid->forces.LE.K;
  const PetscScalar L0  = solid->forces.LE.L0;

  for(PetscInt j=0; j< solid->N_total_lns; j++)      
  {
    const PetscInt p0 = solid->lns[j][0]; 
    const PetscInt p1 = solid->lns[j][1]; 

    CalcNormalVector( solid->x, solid->y, solid->z, p0, p1, &nVx, &nVy, &nVz); 
    L = CalcDistance( solid->x, solid->y, solid->z, p0, p1);

    Ef = K*(L-L0); 

    solid->forces.Fx[p0] += nVx*Ef, solid->forces.Fy[p0] += nVy*Ef, solid->forces.Fz[p0] += nVz*Ef;
    solid->forces.Fx[p1] -= nVx*Ef, solid->forces.Fy[p1] -= nVy*Ef, solid->forces.Fz[p1] -= nVz*Ef;  
  }

  PetscPrintf(PETSC_COMM_WORLD,"  Calculated the linear elastic forces at the solid points based on L & k.\n");
  return 0;
}

PetscErrorCode LinearElasticityOrigin(System_Struct *system, Solid_Struct *solid)
{ 
  PetscScalar nVx=0, nVy=0, nVz=0, L, L0, Ef;
  const PetscScalar K   = solid->forces.LE.K;

  for(PetscInt j=0; j< solid->N_total_lns; j++)      
  {
    const PetscInt p0 = solid->lns[j][0];  
    const PetscInt p1 = solid->lns[j][1];  

    CalcNormalVector( solid->x, solid->y, solid->z, p0, p1, &nVx, &nVy, &nVz); 
    L = CalcDistance( solid->x, solid->y, solid->z, p0, p1);
    L0 = solid->ln_size[j];

    Ef = K*(L-L0); 

    solid->forces.Fx[p0] += nVx*Ef, solid->forces.Fy[p0] += nVy*Ef, solid->forces.Fz[p0] += nVz*Ef;
    solid->forces.Fx[p1] -= nVx*Ef, solid->forces.Fy[p1] -= nVy*Ef, solid->forces.Fz[p1] -= nVz*Ef;  
  }

  PetscPrintf(PETSC_COMM_WORLD,"  Calculated the linear elastic forces at the solid points based on original L0 & K.\n");

  return 0;
}

PetscErrorCode LinearElasticityOriginMultiK(System_Struct *system, Solid_Struct *solid)
{ 
  PetscScalar nVx=0, nVy=0, nVz=0, L, L0, Ef;
  const PetscScalar K   = solid->forces.LE.K;

  for(PetscInt j=0; j< solid->N_total_lns; j++)  
  {
    const PetscScalar M = solid->forces.LE.M[j];    
    const PetscInt p0 = solid->lns[j][0];        
    const PetscInt p1 = solid->lns[j][1];        

    CalcNormalVector( solid->x, solid->y, solid->z, p0, p1, &nVx, &nVy, &nVz); 
    L = CalcDistance( solid->x, solid->y, solid->z, p0, p1);
    L0 = solid->ln_size[j];

    Ef = K*M*(L-L0); 

    
    solid->forces.Fx[p0] += nVx*Ef, solid->forces.Fy[p0] += nVy*Ef, solid->forces.Fz[p0] += nVz*Ef;
    solid->forces.Fx[p1] -= nVx*Ef, solid->forces.Fy[p1] -= nVy*Ef, solid->forces.Fz[p1] -= nVz*Ef;  
  }

  PetscPrintf(PETSC_COMM_WORLD,"  Calculated the linear elastic forces at the solid points based on original L0 & multiple Ks.\n");
  
  return 0;
}

PetscErrorCode LinearElasticityNull(System_Struct *system, Solid_Struct *solid){ 
  return 0;} 

PetscErrorCode LinearElasticityFree(Solid_Struct *solid)
{
  if( solid->forces.LE.Type == FORCE_LE_INITIAL_DISP || 
      solid->forces.LE.Type == FORCE_LE_INITIAL_DISP_MULTI || 
      solid->forces.LE.Type == FORCE_LE_INITIAL_DISP_COND)
      PetscFree(solid->ln_size);

  if(solid->forces.LE.Type == FORCE_LE_INITIAL_DISP_MULTI)
    PetscFree(solid->forces.LE.M);

  return 0;
}

PetscErrorCode LinearElasticityInitialize(Solid_Struct *solid)
{

  if(solid->forces.LE.Type == FORCE_LE_INITIAL_DISP_MULTI)  
    if(solid->forces.LE.M_size != solid->N_total_lns) 
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong size of the Multiplier M"); 

  if (solid->forces.LE.Type == FORCE_LE_INITIAL_DISP_COND)
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! FORCE_LE_INITIAL_DISP_COND not defined"); 

  if( solid->forces.LE.Type == FORCE_LE_INITIAL_DISP || 
      solid->forces.LE.Type == FORCE_LE_INITIAL_DISP_MULTI )  
    SolidCalculateLineSize(solid);  

  return 0; 
}

#endif