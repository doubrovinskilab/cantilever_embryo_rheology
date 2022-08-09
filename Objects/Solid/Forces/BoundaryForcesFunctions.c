#ifndef Solid_Forces_BoundaryForce_Functions_C
#define Solid_Forces_BoundaryForce_Functions_C


PetscErrorCode BoundaryForce_DistFixed(Solid_Struct *solid)
{
  
  for (PetscInt bp = 0; bp<solid->N_bound_pts; bp++)
  {
    const PetscInt g = solid->bound_group_pts[bp];  
    if (solid->Gr_type[g] == DIST_FIXED)
    {
      const PetscInt p = solid->bound_pts[bp];         
      solid->forces.Fx[p] = 0; 
      solid->forces.Fy[p] = 0; 
      solid->forces.Fz[p] = 0; 
    }
  }

  
  for (PetscInt bl = 0; bl<solid->N_bound_lns; bl++)
  {
    const PetscInt g = solid->bound_group_lns[bl];  
    if (solid->Gr_type[g] == DIST_FIXED)
    {
      const PetscInt l = solid->bound_lns[bl];   
      for (PetscInt i=0; i<2; i++){
        solid->forces.Fx[solid->lns[l][i]] = 0;  
        solid->forces.Fy[solid->lns[l][i]] = 0;  
        solid->forces.Fz[solid->lns[l][i]] = 0; }
    }
  }

  
  for (PetscInt bt = 0; bt<solid->N_bound_trs; bt++)
  {
    const PetscInt g = solid->bound_group_trs[bt];  
    if (solid->Gr_type[g] == DIST_FIXED)
    {
      const PetscInt t = solid->bound_trs[bt];   
      for (PetscInt i=0; i<3; i++){
        solid->forces.Fx[solid->trs[t][i]] = 0;  
        solid->forces.Fy[solid->trs[t][i]] = 0;  
        solid->forces.Fz[solid->trs[t][i]] = 0; }
    }
  }

  return 0;
}

PetscErrorCode BoundaryForce_DistFixedNull(Solid_Struct *solid) { return 0;}


PetscErrorCode BoundaryForce_InsertForce(Solid_Struct *solid)
{
  PetscScalar *Gr_size;                 
  PetscMalloc1(solid->N_Gr, &Gr_size);
  for(PetscInt g=0; g<solid->N_Gr; g++)
  {
    Gr_size[g] = 0.0;
    if(solid->Gr_elem_type[g]==POINT)
    {
      Gr_size[g] = 1.0;
    }
    else if(solid->Gr_elem_type[g]==LINE) 
    {
      for(PetscInt bl=0; bl<solid->N_bound_lns; bl++) 
      {
        if(solid->bound_group_lns[bl] != g)   
          continue;
        const PetscInt l = solid->bound_group_lns[bl];
        Gr_size[g] += CalcDistance(solid->x,solid->y,solid->z,solid->lns[l][0],solid->lns[l][1]);
      }
    }
    else if(solid->Gr_elem_type[g]==TRIANGLE) 
    {
      for(PetscInt bt=0; bt<solid->N_bound_trs; bt++) 
      {
        if(solid->bound_group_trs[bt] != g)  
          continue;
        const PetscInt t = solid->bound_group_trs[bt];
        Gr_size[g] += CalcAreaOfTriangle(solid->x,solid->y,solid->z,solid->trs[t][0],solid->trs[t][1],solid->trs[t][2]);
      }
    }
  }    

  for (PetscInt bp = 0; bp<solid->N_bound_pts; bp++)
  {
    const PetscInt g = solid->bound_group_pts[bp];  
    if (solid->Gr_type[g] == FORCE_INSERT)
    {
      const PetscInt p = solid->bound_pts[bp];         
      const PetscScalar BFx = solid->Gr_value[g][0];   
      const PetscScalar BFy = solid->Gr_value[g][1];   
      const PetscScalar BFz = solid->Gr_value[g][2];   

      solid->forces.Fx[p] = BFx; 
      solid->forces.Fy[p] = BFy; 
      solid->forces.Fz[p] = BFz; 
    }
  }
  
  for (PetscInt bl = 0; bl<solid->N_bound_lns; bl++)
  {
    PetscScalar bfx=0, bfy=0, bfz=0, bs=0;   
    const PetscInt g = solid->bound_group_lns[bl];  
    if (solid->Gr_type[g] == FORCE_INSERT)
    {
      const PetscInt l = solid->bound_lns[bl];   
      const PetscScalar D   = CalcDistance(solid->x, solid->y, solid->z, solid->lns[l][0], solid->lns[l][1]); 
      const PetscScalar BFx = solid->Gr_value[g][0];   
      const PetscScalar BFy = solid->Gr_value[g][1];   
      const PetscScalar BFz = solid->Gr_value[g][2];   
      const PetscScalar BS  = Gr_size[g];              

      bfx=BFx/BS, bfy=BFy/BS, bfz=BFz/BS;   
      bs = D/2.0;                           

      for (PetscInt i=0; i<2; i++){
        solid->forces.Fx[solid->lns[l][i]] = bfx*bs;  
        solid->forces.Fy[solid->lns[l][i]] = bfy*bs;  
        solid->forces.Fz[solid->lns[l][i]] = bfz*bs; }
    }
  }

  for (PetscInt bt = 0; bt<solid->N_bound_trs; bt++)
  {
    PetscScalar bfx=0, bfy=0, bfz=0, bs=0;   
    const PetscInt g = solid->bound_group_trs[bt];  
    if (solid->Gr_type[g] == FORCE_INSERT)
    {
      const PetscInt t = solid->bound_trs[bt];   
      const PetscScalar A   = CalcAreaOfTriangle(solid->x, solid->y, solid->z, solid->trs[t][0], solid->trs[t][1], solid->trs[t][2]); 
      const PetscScalar BFx = solid->Gr_value[g][0];   
      const PetscScalar BFy = solid->Gr_value[g][1];   
      const PetscScalar BFz = solid->Gr_value[g][2];   
      const PetscScalar BS  = Gr_size[g];              

      bfx=BFx/BS, bfy=BFy/BS, bfz=BFz/BS;   
      bs = A/3.0;                           

      for (PetscInt i=0; i<3; i++){
        solid->forces.Fx[solid->trs[t][i]] = bfx*bs;  
        solid->forces.Fy[solid->trs[t][i]] = bfy*bs;  
        solid->forces.Fz[solid->trs[t][i]] = bfz*bs; }
    }
  }
  PetscFree(Gr_size);
  return 0;
}


PetscErrorCode BoundaryForce_InsertForcePoints(Solid_Struct *solid)
{
  for (PetscInt bp = 0; bp<solid->N_bound_pts; bp++)
  {
    const PetscInt g = solid->bound_group_pts[bp];  
    if (solid->Gr_type[g] == FORCE_ADD)
    {
      const PetscInt p = solid->bound_pts[bp];         
      const PetscScalar BFx = solid->Gr_value[g][0];   
      const PetscScalar BFy = solid->Gr_value[g][1];   
      const PetscScalar BFz = solid->Gr_value[g][2];   

      solid->forces.Fx[p] = BFx; 
      solid->forces.Fy[p] = BFy; 
      solid->forces.Fz[p] = BFz; 
    }
  }
  return 0;
}


PetscErrorCode BoundaryForce_InsertForceLines(Solid_Struct *solid)
{
  PetscScalar *Gr_size;                 
  PetscMalloc1(solid->N_Gr, &Gr_size);
  for(PetscInt g=0; g<solid->N_Gr; g++)
  {
    Gr_size[g] = 0.0;
    if(solid->Gr_elem_type[g]==LINE) 
    {
      for(PetscInt bl=0; bl<solid->N_bound_lns; bl++) 
      {
        if(solid->bound_group_lns[bl] != g)   
          continue;
        const PetscInt l = solid->bound_group_lns[bl];
        Gr_size[g] += CalcDistance(solid->x,solid->y,solid->z,solid->lns[l][0],solid->lns[l][1]);
      }
    }
  }    

  for (PetscInt bl = 0; bl<solid->N_bound_lns; bl++)
  {
    PetscScalar bfx=0, bfy=0, bfz=0, bs=0;   
    const PetscInt g = solid->bound_group_lns[bl];  
    if (solid->Gr_type[g] == FORCE_ADD)
    {
      const PetscInt l = solid->bound_lns[bl];          
      const PetscScalar D   = CalcDistance(solid->x, solid->y, solid->z, solid->lns[l][0], solid->lns[l][1]); 
      const PetscScalar BFx = solid->Gr_value[g][0];    
      const PetscScalar BFy = solid->Gr_value[g][1];    
      const PetscScalar BFz = solid->Gr_value[g][2];    
      const PetscScalar BS  = Gr_size[g];               

      bfx=BFx/BS, bfy=BFy/BS, bfz=BFz/BS;   
      bs = D/2.0;                           

      for (PetscInt i=0; i<2; i++){
        solid->forces.Fx[solid->lns[l][i]] = bfx*bs;  
        solid->forces.Fy[solid->lns[l][i]] = bfy*bs;  
        solid->forces.Fz[solid->lns[l][i]] = bfz*bs; }
    }
  }  
  PetscFree(Gr_size);
  return 0;
}

PetscErrorCode BoundaryForce_InsertForceTriangles(Solid_Struct *solid)
{
  PetscScalar *Gr_size;                 
  PetscMalloc1(solid->N_Gr, &Gr_size);
  for(PetscInt g=0; g<solid->N_Gr; g++)
  {
    if(solid->Gr_elem_type[g]==TRIANGLE) 
    {
      for(PetscInt bt=0; bt<solid->N_bound_trs; bt++) 
      {
        if(solid->bound_group_trs[bt] != g)  
          continue;
        const PetscInt t = solid->bound_group_trs[bt];
        Gr_size[g] += CalcAreaOfTriangle(solid->x,solid->y,solid->z,solid->trs[t][0],solid->trs[t][1],solid->trs[t][2]);
      }
    }
  }
  
  for (PetscInt bt = 0; bt<solid->N_bound_trs; bt++)
  {
    PetscScalar bfx=0, bfy=0, bfz=0, bs=0;   
    const PetscInt g = solid->bound_group_trs[bt];  
    if (solid->Gr_type[g] == FORCE_ADD)
    {
      const PetscInt t = solid->bound_trs[bt];          
      const PetscScalar A   = CalcAreaOfTriangle(solid->x, solid->y, solid->z, solid->trs[t][0], solid->trs[t][1], solid->trs[t][2]); 
      const PetscScalar BFx = solid->Gr_value[g][0];    
      const PetscScalar BFy = solid->Gr_value[g][1];    
      const PetscScalar BFz = solid->Gr_value[g][2];    
      const PetscScalar BS  = Gr_size[g];               

      bfx=BFx/BS, bfy=BFy/BS, bfz=BFz/BS;   
      bs = A/3.0;                           

      for (PetscInt i=0; i<3; i++){
        solid->forces.Fx[solid->trs[t][i]] = bfx*bs;  
        solid->forces.Fy[solid->trs[t][i]] = bfy*bs;  
        solid->forces.Fz[solid->trs[t][i]] = bfz*bs; }
    }
  }
  PetscFree(Gr_size);
  return 0;
}

PetscErrorCode BoundaryForce_InsertForceNull(Solid_Struct *solid) {return 0;}


PetscErrorCode BoundaryForce_AddForce(Solid_Struct *solid)
{
  PetscScalar *Gr_size;                 
  PetscMalloc1(solid->N_Gr, &Gr_size);
  for(PetscInt g=0; g<solid->N_Gr; g++)
  {
    Gr_size[g] = 0.0;
    if(solid->Gr_elem_type[g]==POINT)
    {
      Gr_size[g] = 1.0;
    }
    else if(solid->Gr_elem_type[g]==LINE) 
    {
      for(PetscInt bl=0; bl<solid->N_bound_lns; bl++) 
      {
        if(solid->bound_group_lns[bl] != g)   
          continue;
        const PetscInt l = solid->bound_group_lns[bl];
        Gr_size[g] += CalcDistance(solid->x,solid->y,solid->z,solid->lns[l][0],solid->lns[l][1]);
      }
    }
    else if(solid->Gr_elem_type[g]==TRIANGLE) 
    {
      for(PetscInt bt=0; bt<solid->N_bound_trs; bt++) 
      {
        if(solid->bound_group_trs[bt] != g)  
          continue;
        const PetscInt t = solid->bound_group_trs[bt];
        Gr_size[g] += CalcAreaOfTriangle(solid->x,solid->y,solid->z,solid->trs[t][0],solid->trs[t][1],solid->trs[t][2]);
      }
    }
  }    

  
  for (PetscInt bp = 0; bp<solid->N_bound_pts; bp++)
  {
    const PetscInt g = solid->bound_group_pts[bp];  
    if (solid->Gr_type[g] == FORCE_ADD)
    {
      const PetscInt p = solid->bound_pts[bp];         
      const PetscScalar BFx = solid->Gr_value[g][0];   
      const PetscScalar BFy = solid->Gr_value[g][1];   
      const PetscScalar BFz = solid->Gr_value[g][2];   

      solid->forces.Fx[p] += BFx; 
      solid->forces.Fy[p] += BFy; 
      solid->forces.Fz[p] += BFz; 
    }
  }

  
  for (PetscInt bl = 0; bl<solid->N_bound_lns; bl++)
  {
    PetscScalar bfx=0, bfy=0, bfz=0, bs=0;   
    const PetscInt g = solid->bound_group_lns[bl];  
    if (solid->Gr_type[g] == FORCE_ADD)
    {
      const PetscInt l = solid->bound_lns[bl];          
      const PetscScalar D   = CalcDistance(solid->x, solid->y, solid->z, solid->lns[l][0], solid->lns[l][1]); 
      const PetscScalar BFx = solid->Gr_value[g][0];    
      const PetscScalar BFy = solid->Gr_value[g][1];    
      const PetscScalar BFz = solid->Gr_value[g][2];    
      const PetscScalar BS  = Gr_size[g];               

      bfx=BFx/BS, bfy=BFy/BS, bfz=BFz/BS;   
      bs = D/2.0;                           

      for (PetscInt i=0; i<2; i++){
        solid->forces.Fx[solid->lns[l][i]] += bfx*bs;  
        solid->forces.Fy[solid->lns[l][i]] += bfy*bs;  
        solid->forces.Fz[solid->lns[l][i]] += bfz*bs; }
    }
  }

  
  for (PetscInt bt = 0; bt<solid->N_bound_trs; bt++)
  {
    PetscScalar bfx=0, bfy=0, bfz=0, bs=0;   
    const PetscInt g = solid->bound_group_trs[bt];  
    if (solid->Gr_type[g] == FORCE_ADD)
    {
      const PetscInt t = solid->bound_trs[bt];          
      const PetscScalar A   = CalcAreaOfTriangle(solid->x, solid->y, solid->z, solid->trs[t][0], solid->trs[t][1], solid->trs[t][2]); 
      const PetscScalar BFx = solid->Gr_value[g][0];    
      const PetscScalar BFy = solid->Gr_value[g][1];    
      const PetscScalar BFz = solid->Gr_value[g][2];    
      const PetscScalar BS  = Gr_size[g];               

      bfx=BFx/BS, bfy=BFy/BS, bfz=BFz/BS;   
      bs = A/3.0;                           

      for (PetscInt i=0; i<3; i++){
        solid->forces.Fx[solid->trs[t][i]] += bfx*bs;  
        solid->forces.Fy[solid->trs[t][i]] += bfy*bs;  
        solid->forces.Fz[solid->trs[t][i]] += bfz*bs; }
    }
  }
  PetscFree(Gr_size);
  return 0;
}


PetscErrorCode BoundaryForce_AddForcePoints(Solid_Struct *solid)
{

  for (PetscInt bp = 0; bp<solid->N_bound_pts; bp++)
  {
    const PetscInt g = solid->bound_group_pts[bp];  
    if (solid->Gr_type[g] == FORCE_ADD)
    {
      const PetscInt p = solid->bound_pts[bp];         
      const PetscScalar BFx = solid->Gr_value[g][0];   
      const PetscScalar BFy = solid->Gr_value[g][1];   
      const PetscScalar BFz = solid->Gr_value[g][2];   

      solid->forces.Fx[p] += BFx; 
      solid->forces.Fy[p] += BFy; 
      solid->forces.Fz[p] += BFz; 
    }
  }
  return 0;
}


PetscErrorCode BoundaryForce_AddForceLines(Solid_Struct *solid)
{
  PetscScalar *Gr_size;                 
  
  PetscMalloc1(solid->N_Gr, &Gr_size);
  for(PetscInt g=0; g<solid->N_Gr; g++)
  {
    Gr_size[g] = 0.0;
    if(solid->Gr_elem_type[g]==LINE) 
    {
      for(PetscInt bl=0; bl<solid->N_bound_lns; bl++) 
      {
        if(solid->bound_group_lns[bl] != g)   
          continue;
        const PetscInt l = solid->bound_group_lns[bl];
        Gr_size[g] += CalcDistance(solid->x,solid->y,solid->z,solid->lns[l][0],solid->lns[l][1]);
      }
    }
  }    

  
  for (PetscInt bl = 0; bl<solid->N_bound_lns; bl++)
  {
    PetscScalar bfx=0, bfy=0, bfz=0, bs=0;   
    const PetscInt g = solid->bound_group_lns[bl];  
    if (solid->Gr_type[g] == FORCE_ADD)
    {
      const PetscInt l = solid->bound_lns[bl];          
      const PetscScalar D   = CalcDistance(solid->x, solid->y, solid->z, solid->lns[l][0], solid->lns[l][1]); 
      const PetscScalar BFx = solid->Gr_value[g][0];    
      const PetscScalar BFy = solid->Gr_value[g][1];    
      const PetscScalar BFz = solid->Gr_value[g][2];    
      const PetscScalar BS  = Gr_size[g];               

      bfx=BFx/BS, bfy=BFy/BS, bfz=BFz/BS;   
      bs = D/2.0;                           

      for (PetscInt i=0; i<2; i++){
        solid->forces.Fx[solid->lns[l][i]] += bfx*bs;  
        solid->forces.Fy[solid->lns[l][i]] += bfy*bs;  
        solid->forces.Fz[solid->lns[l][i]] += bfz*bs; }
    }
  }
  PetscFree(Gr_size);
  return 0;
}

PetscErrorCode BoundaryForce_AddForceTriangles(Solid_Struct *solid)
{
  PetscScalar *Gr_size;                 

  PetscMalloc1(solid->N_Gr, &Gr_size);
  for(PetscInt g=0; g<solid->N_Gr; g++)
  {
    if(solid->Gr_elem_type[g]==TRIANGLE) 
    {
      for(PetscInt bt=0; bt<solid->N_bound_trs; bt++) 
      {
        if(solid->bound_group_trs[bt] != g)  
          continue;
        const PetscInt t = solid->bound_group_trs[bt];
        Gr_size[g] += CalcAreaOfTriangle(solid->x,solid->y,solid->z,solid->trs[t][0],solid->trs[t][1],solid->trs[t][2]);
      }
    }
  }

  
  for (PetscInt bt = 0; bt<solid->N_bound_trs; bt++)
  {
    PetscScalar bfx=0, bfy=0, bfz=0, bs=0;   
    const PetscInt g = solid->bound_group_trs[bt];  
    if (solid->Gr_type[g] == FORCE_ADD)
    {
      const PetscInt t = solid->bound_trs[bt];          
      const PetscScalar A   = CalcAreaOfTriangle(solid->x, solid->y, solid->z, solid->trs[t][0], solid->trs[t][1], solid->trs[t][2]); 
      const PetscScalar BFx = solid->Gr_value[g][0];    
      const PetscScalar BFy = solid->Gr_value[g][1];    
      const PetscScalar BFz = solid->Gr_value[g][2];    
      const PetscScalar BS  = Gr_size[g];               

      bfx=BFx/BS, bfy=BFy/BS, bfz=BFz/BS;   
      bs = A/3.0;                           

      for (PetscInt i=0; i<3; i++){
        solid->forces.Fx[solid->trs[t][i]] += bfx*bs;  
        solid->forces.Fy[solid->trs[t][i]] += bfy*bs;  
        solid->forces.Fz[solid->trs[t][i]] += bfz*bs; }
    }
  }
  PetscFree(Gr_size);
  return 0;
}

PetscErrorCode BoundaryForce_AddForceNull(Solid_Struct *solid) {return 0;}


PetscErrorCode BoundaryForce_AddForceEllipsePoints(Solid_Struct *solid)
{
  
  PetscScalar Theta, nVx, nVy, nVz;
  for (PetscInt bp = 0; bp<solid->N_bound_pts; bp++)
  {
    const PetscInt g = solid->bound_group_pts[bp];  
    if (solid->Gr_type[g] == FORCE_ADD_ELLIPSE)
    {
      const PetscInt p = solid->bound_pts[bp]; 

      const PetscScalar A   = solid->Gr_value[g][0];   
      const PetscScalar B   = solid->Gr_value[g][1];   
      const PetscInt Dir    = (PetscInt) solid->Gr_value[g][2];  
      const PetscScalar Mag = solid->Gr_value[g][3];  

      Theta = CalcEllipsoidAngle(solid->x[p], solid->y[p], A, B);

      CalcNormalVector2( solid->x[p], solid->y[p], 0, A*PetscCosScalar(Theta+Dir*1e-6), B*PetscSinScalar(Theta+Dir*1e-6), 0, &nVx, &nVy, &nVz);

      solid->forces.Fx[p] += Mag*nVx; 
      solid->forces.Fy[p] += Mag*nVy; 
      solid->forces.Fz[p] += 0; 
    }
  }
  return 0;
}


PetscErrorCode BoundaryForce_AddForceEllipseNull(Solid_Struct *solid) {return 0;}

PetscErrorCode BoundaryForce_SpeedForceReload(Solid_Struct *solid, const PetscScalar reload_time)
{
  if(solid->forces.BC_FS.N>0)
  { 
    for (PetscInt i=0; i< solid->forces.BC_FS.N; i++)         
    {
      const PetscInt g = solid->forces.BC_FS.g[i];

      
      solid->forces.BC_FS.x[i] += solid->Gr_value[g][0]*reload_time;
      solid->forces.BC_FS.y[i] += solid->Gr_value[g][1]*reload_time;
      solid->forces.BC_FS.z[i] += solid->Gr_value[g][2]*reload_time;
    }
  }
  return 0;
}

PetscErrorCode BoundaryForce_SpeedForceInitialize(Solid_Struct *solid)
{

  PetscBool Found;
  
  solid->forces.BC_FS.N = 0;                        

  for(PetscInt j=0; j< solid->N_bound_pts; j++) {   
    const PetscInt g = solid->bound_group_pts[j];   
    if (solid->Gr_type[g] == FORCE_SPEED) solid->forces.BC_FS.N++; }

  for(PetscInt j=0; j< solid->N_bound_lns; j++) {   
    const PetscInt g = solid->bound_group_lns[j];   
    if (solid->Gr_type[g] == FORCE_SPEED) solid->forces.BC_FS.N+=2; }

  for(PetscInt j=0; j< solid->N_bound_trs; j++) {   
    const PetscInt g = solid->bound_group_trs[j];   
    if (solid->Gr_type[g] == FORCE_SPEED) solid->forces.BC_FS.N+=3; }
  
  if (solid->forces.BC_FS.N > 0) 
  {
    PetscMalloc5(solid->forces.BC_FS.N, &solid->forces.BC_FS.p, solid->forces.BC_FS.N, &solid->forces.BC_FS.g, 
      solid->forces.BC_FS.N, &solid->forces.BC_FS.x, solid->forces.BC_FS.N, &solid->forces.BC_FS.y, solid->forces.BC_FS.N, &solid->forces.BC_FS.z);

    for(PetscInt j=0; j< solid->forces.BC_FS.N; j++) solid->forces.BC_FS.p[j] = -1;

    PetscInt k = 0;
    for(PetscInt j=0; j< solid->N_bound_pts; j++) 
    {
      const PetscInt g = solid->bound_group_pts[j]; 
      if (solid->Gr_type[g] == FORCE_SPEED) {
        solid->forces.BC_FS.p[k] = solid->bound_pts[j];
        solid->forces.BC_FS.g[k] = solid->bound_group_pts[j];
        k++;}
    }

    for(PetscInt j=0; j< solid->N_bound_lns; j++) 
    {
      const PetscInt g = solid->bound_group_lns[j]; 
      const PetscInt l = solid->bound_lns[j];       
      if (solid->Gr_type[g] == FORCE_SPEED) 
        for (PetscInt i=0; i<2; i++) 
        {
          const PetscInt p = solid->lns[l][i]; 
          Found = PETSC_FALSE;
          for (PetscInt q = 0; q < k ; q++)
            if ( solid->forces.BC_FS.p[q] == p){
              Found = PETSC_TRUE;
              break;
            }
          
          if (PetscNot(Found)){
            solid->forces.BC_FS.p[k] = p;
            solid->forces.BC_FS.g[k] = g;
            k++;}
        }
    }

    for(PetscInt j=0; j< solid->N_bound_trs; j++) 
    {  
      const PetscInt g = solid->bound_group_trs[j]; 
      const PetscInt t = solid->bound_trs[j]      ; 
      if (solid->Gr_type[g] == FORCE_SPEED) 
        for (PetscInt i=0; i<3; i++) 
        {
          const PetscInt p = solid->trs[t][i]; 
          Found = PETSC_FALSE;
          for (PetscInt q = 0; q < k ; q++)
            if ( solid->forces.BC_FS.p[q] == p){
              Found = PETSC_TRUE;
              break;
            }
          
          if (PetscNot(Found)){
            solid->forces.BC_FS.p[k] = p;
            solid->forces.BC_FS.g[k] = g;
            k++;}
        }
    }

    if(k>solid->forces.BC_FS.N)
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error in create_points_for_BC_FORCE_SPEED");  
    else
      solid->forces.BC_FS.N = k;
    

    for(PetscInt j=0; j< solid->forces.BC_FS.N; j++)
    {
      const PetscInt p = solid->forces.BC_FS.p[j];
      const PetscInt g = solid->forces.BC_FS.g[j]; 
      solid->forces.BC_FS.x[j] = solid->x[p] + solid->Gr_value[g][3];
      solid->forces.BC_FS.y[j] = solid->y[p] + solid->Gr_value[g][4];
      solid->forces.BC_FS.z[j] = solid->z[p] + solid->Gr_value[g][5];
    }
  }

  return 0;
}

PetscErrorCode BoundaryForce_SpeedForce(System_Struct *system, Solid_Struct *solid)
{
  PetscScalar External_Force_x = 0;
  PetscScalar External_Force_y = 0;
  PetscScalar External_Force_z = 0;
  PetscViewer  viewer; 
  PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERASCII);
  PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);
  PetscViewerFileSetName(viewer, solid->forces.BC_FS.File); 

  for (PetscInt i=0; i< solid->forces.BC_FS.N; i++)         
  {
    const PetscInt p = solid->forces.BC_FS.p[i];  
    const PetscInt g = solid->forces.BC_FS.g[i];  
    const PetscScalar K  = solid->Gr_value[g][6];    
    const PetscScalar Vx = solid->forces.BC_FS.x[i] - solid->x[p]; 
    const PetscScalar Vy = solid->forces.BC_FS.y[i] - solid->y[p]; 
    const PetscScalar Vz = solid->forces.BC_FS.z[i] - solid->z[p];     
    External_Force_x += -1*solid->forces.Fx[p] + K*Vx;
    External_Force_y += -1*solid->forces.Fy[p] + K*Vy;
    External_Force_z += -1*solid->forces.Fz[p] + K*Vz;    
    solid->forces.Fx[p] = K*Vx; 
    solid->forces.Fy[p] = K*Vy; 
    solid->forces.Fz[p] = K*Vz; 
  }
  PetscViewerASCIIPrintf(viewer, "%g %g %g %g \n", (double)system->time, (double)External_Force_x, (double)External_Force_y, (double)External_Force_z);
  PetscViewerDestroy(&viewer);
  return 0;
}

PetscErrorCode BoundaryForce_SpeedForceNull(System_Struct *system, Solid_Struct *solid) { return 0;}

PetscErrorCode BoundaryForce_SpeedForceEllipseReload(Solid_Struct *solid, const PetscScalar reload_time)
{
  if(solid->forces.BC_FSE.N>0)
  { 
    for(PetscInt j=0; j< solid->forces.BC_FSE.N; j++)
    {
      const PetscInt g = solid->forces.BC_FSE.g[j];
      const PetscScalar DThetaDt  = solid->Gr_value[g][2];
      solid->forces.BC_FSE.theta[j] += DThetaDt*reload_time; 
    }
  }
  return 0;
}

PetscErrorCode BoundaryForce_SpeedForceEllipseInitialize(Solid_Struct *solid)
{
  
  PetscBool Found;
  solid->forces.BC_FSE.N = 0;                       

  for(PetscInt j=0; j< solid->N_bound_pts; j++) {   
    const PetscInt g = solid->bound_group_pts[j];   
    if (solid->Gr_type[g] == FORCE_SPEED_ELLIPSE)   solid->forces.BC_FSE.N++; }

  for(PetscInt j=0; j< solid->N_bound_lns; j++) {   
    const PetscInt g = solid->bound_group_lns[j];   
    if (solid->Gr_type[g] == FORCE_SPEED_ELLIPSE)   solid->forces.BC_FSE.N+=2; }

  for(PetscInt j=0; j< solid->N_bound_trs; j++) {   
    const PetscInt g = solid->bound_group_trs[j];   
    if (solid->Gr_type[g] == FORCE_SPEED_ELLIPSE)   solid->forces.BC_FSE.N+=3; }

  if (solid->forces.BC_FSE.N > 0) 
  {
    PetscMalloc3( solid->forces.BC_FSE.N, &solid->forces.BC_FSE.p, 
                  solid->forces.BC_FSE.N, &solid->forces.BC_FSE.g, 
                  solid->forces.BC_FSE.N, &solid->forces.BC_FSE.theta);
    
    for(PetscInt j=0; j< solid->forces.BC_FSE.N ; j++) 
      solid->forces.BC_FSE.p[j] = -1, solid->forces.BC_FSE.g[j] = -1, solid->forces.BC_FSE.theta[j] =0;

    PetscInt k = 0;
    for(PetscInt j=0; j< solid->N_bound_pts; j++) 
    {
      const PetscInt g = solid->bound_group_pts[j]; 
      if (solid->Gr_type[g] == FORCE_SPEED_ELLIPSE) {
        solid->forces.BC_FSE.p[k] = solid->bound_pts[j];
        solid->forces.BC_FSE.g[k] = solid->bound_group_pts[j];
        k++;}
    }

    for(PetscInt j=0; j< solid->N_bound_lns; j++) 
    {
      const PetscInt g = solid->bound_group_lns[j]; 
      const PetscInt l = solid->bound_lns[j];       
      if (solid->Gr_type[g] == FORCE_SPEED_ELLIPSE) 
        for (PetscInt i=0; i<2; i++) 
        {
          const PetscInt p = solid->lns[l][i]; 
          Found = PETSC_FALSE;
          for (PetscInt q = 0; q < k ; q++)
            if ( solid->forces.BC_FSE.p[q] == p){
              Found = PETSC_TRUE;
              break;
            }
          
          if (PetscNot(Found)){
            solid->forces.BC_FSE.p[k] = p;
            solid->forces.BC_FSE.g[k] = g;
            k++;}
        }
    }

    for(PetscInt j=0; j< solid->N_bound_trs; j++) 
    {  
      const PetscInt g = solid->bound_group_trs[j]; 
      const PetscInt t = solid->bound_trs[j]      ; 
      if (solid->Gr_type[g] == FORCE_SPEED_ELLIPSE) 
        for (PetscInt i=0; i<3; i++) 
        {
          const PetscInt p = solid->trs[t][i]; 
          Found = PETSC_FALSE;
          for (PetscInt q = 0; q < k ; q++)
            if ( solid->forces.BC_FSE.p[q] == p){
              Found = PETSC_TRUE;
              break;
            }
          
          if (PetscNot(Found)){
            solid->forces.BC_FSE.p[k] = p;
            solid->forces.BC_FSE.g[k] = g;
            k++;}
        }
    }

    if(k>solid->forces.BC_FSE.N)
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error in create_points_for_BC_FORCE_SPEED_ELLIPSE");  
    else
      solid->forces.BC_FSE.N = k; 
    
    for(PetscInt j=0; j< solid->forces.BC_FSE.N; j++)
    {
      const PetscInt p = solid->forces.BC_FSE.p[j];
      const PetscInt g = solid->forces.BC_FSE.g[j];
      const PetscScalar A       = solid->Gr_value[g][0];
      const PetscScalar B       = solid->Gr_value[g][1];
      const PetscScalar DTheta0 = solid->Gr_value[g][5];
      const PetscScalar x = solid->x[p];
      const PetscScalar y = solid->y[p];
      PetscScalar Theta = CalcEllipsoidAngle(x, y, A, B);
      solid->forces.BC_FSE.theta[j] = Theta + DTheta0;
    }
  }
  return 0;
}


PetscErrorCode BoundaryForce_SpeedForceEllipse(System_Struct *system, Solid_Struct *solid)
{
  PetscScalar External_Force_x = 0;
  PetscScalar External_Force_y = 0;
  PetscScalar External_Force_z = 0;

  PetscViewer  viewer; 

  PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERASCII);
  PetscViewerFileSetMode(viewer, FILE_MODE_APPEND);
  PetscViewerFileSetName(viewer, solid->forces.BC_FSE.File); 

  for (PetscInt i=0; i< solid->forces.BC_FSE.N; i++)         
  {
    const PetscInt p          = solid->forces.BC_FSE.p[i];      
    const PetscInt g          = solid->forces.BC_FSE.g[i];      
    const PetscScalar Theta   = solid->forces.BC_FSE.theta[i];  
    const PetscScalar A         = solid->Gr_value[g][0];
    const PetscScalar B         = solid->Gr_value[g][1];
    const PetscScalar K         = solid->Gr_value[g][3];
    const PetscScalar Z_Plane   = solid->Gr_value[g][4];    
    const PetscScalar Vx = A*PetscCosScalar(Theta) - solid->x[p];
    const PetscScalar Vy = B*PetscSinScalar(Theta) - solid->y[p];
    const PetscScalar Vz = Z_Plane - solid->z[p];    
    External_Force_x += -1*solid->forces.Fx[p] + K*Vx; 
    External_Force_y += -1*solid->forces.Fy[p] + K*Vy; 
    External_Force_z += -1*solid->forces.Fz[p] + K*Vz; 
    solid->forces.Fx[p] = K*Vx; 
    solid->forces.Fy[p] = K*Vy; 
    solid->forces.Fz[p] = K*Vz; 
  }
  PetscViewerASCIIPrintf(viewer, "%g %g %g %g \n", (double)system->time, (double)External_Force_x, (double)External_Force_y, (double)External_Force_z);
  PetscViewerDestroy(&viewer);
  return 0;
}


PetscErrorCode BoundaryForce_SpeedForceEllipseNull(System_Struct *system, Solid_Struct *solid) { return 0;}

#endif 
