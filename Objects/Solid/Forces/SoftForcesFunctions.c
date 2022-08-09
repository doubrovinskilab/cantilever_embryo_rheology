#ifndef Solid_Forces_SoftForce_Functions_C
#define Solid_Forces_SoftForce_Functions_C

PetscErrorCode SoftForceCreateOneEllipsoid(System_Struct *system, Solid_Struct *solid)
{
  const PetscInt Dim = system->Dim;
  const PetscScalar ds = solid->forces.SF.ds;
  const PetscScalar As = solid->forces.SF.A - ds;
  const PetscScalar Al = solid->forces.SF.A + ds;
  const PetscScalar Bs = solid->forces.SF.B - ds;
  const PetscScalar Bl = solid->forces.SF.B + ds;
  const PetscScalar Cs = solid->forces.SF.C - ds;
  const PetscScalar Cl = solid->forces.SF.C + ds;
  PetscScalar Es, El = 0; 
  solid->forces.SF.N = 0;
  for (PetscInt p=0; p<solid->N_total_pts; p++){
    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    if (Dim == 3) 
      Es = x*x/(As*As) + y*y/(Bs*Bs) + z*z/(Cs*Cs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl) + z*z/(Cl*Cl);
    else          
      Es = x*x/(As*As) + y*y/(Bs*Bs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl);

    if (Es>=1.0 && El<=1.0) solid->forces.SF.N++;
  }

  PetscMalloc1(solid->forces.SF.N, &solid->forces.SF.p);
  PetscInt j = 0;
  for (PetscInt p=0; p<solid->N_total_pts; p++){
    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    if (Dim == 3) 
      Es = x*x/(As*As) + y*y/(Bs*Bs) + z*z/(Cs*Cs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl) + z*z/(Cl*Cl);
    else          
      Es = x*x/(As*As) + y*y/(Bs*Bs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl);

    if (Es>=1.0 && El<=1.0) {
      solid->forces.SF.p[j] = p;
      j++;}
  }

  PetscScalar Rho;

  const PetscScalar A  = solid->forces.SF.A;
  const PetscScalar B  = solid->forces.SF.B;
  const PetscScalar C  = solid->forces.SF.C;

  for (PetscInt j=0; j<solid->forces.SF.N; j++)
  {
    const PetscInt p = solid->forces.SF.p[j];

    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    Rho = PetscSqrtScalar( (x*x)/(A*A) + (y*y)/(B*B) + (z*z)/(C*C) );

    solid->x[p] = x/Rho;
    solid->y[p] = y/Rho;
    solid->z[p] = z/Rho*(Dim-2);
  }
  return 0;
}


PetscErrorCode SoftForceCreateTwoEllipsoids(System_Struct *system, Solid_Struct *solid)
{
  const PetscInt Dim = system->Dim;
  const PetscScalar ds  = solid->forces.SF.ds;
  const PetscScalar As  = solid->forces.SF.A  - ds;
  const PetscScalar Al  = solid->forces.SF.A  + ds;
  const PetscScalar Bds = solid->forces.SF.Bd - ds;
  const PetscScalar Bdl = solid->forces.SF.Bd + ds;
  const PetscScalar Bvs = solid->forces.SF.Bv - ds;
  const PetscScalar Bvl = solid->forces.SF.Bv + ds;
  const PetscScalar Cs  = solid->forces.SF.C  - ds;
  const PetscScalar Cl  = solid->forces.SF.C  + ds;
  PetscScalar Bs, Bl;
  PetscScalar Es, El = 0; 
 
  solid->forces.SF.N = 0;
  for (PetscInt p=0; p<solid->N_total_pts; p++){
    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    if (y>0) Bs = Bds, Bl = Bdl;
    else     Bs = Bvs, Bl = Bvl;

    if (Dim == 3) 
      Es = x*x/(As*As) + y*y/(Bs*Bs) + z*z/(Cs*Cs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl) + z*z/(Cl*Cl);
    else          
      Es = x*x/(As*As) + y*y/(Bs*Bs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl);

    if (Es>=1.0 && El<=1.0) solid->forces.SF.N++;
  }

  PetscMalloc1(solid->forces.SF.N, &solid->forces.SF.p);

  
  PetscInt j = 0;
  for (PetscInt p=0; p<solid->N_total_pts; p++){
    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    if (y>0) Bs = Bds, Bl = Bdl;
    else     Bs = Bvs, Bl = Bvl;

    if (Dim == 3) 
      Es = x*x/(As*As) + y*y/(Bs*Bs) + z*z/(Cs*Cs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl) + z*z/(Cl*Cl);
    else          
      Es = x*x/(As*As) + y*y/(Bs*Bs),
      El = x*x/(Al*Al) + y*y/(Bl*Bl);

    if (Es>=1.0 && El<=1.0) {
      solid->forces.SF.p[j] = p;
      j++;}
  }

  PetscScalar Rho;

  const PetscScalar A  = solid->forces.SF.A;
  const PetscScalar Bd = solid->forces.SF.Bd;
  const PetscScalar Bv = solid->forces.SF.Bv;
  const PetscScalar C  = solid->forces.SF.C;
  PetscScalar B;

  for (PetscInt j=0; j<solid->forces.SF.N; j++)
  {
    const PetscInt p = solid->forces.SF.p[j];

    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    if (y>0) B = Bd;
    else     B = Bv;

    Rho = PetscSqrtScalar( (x*x)/(A*A) + (y*y)/(B*B) + (z*z)/(C*C) );

    solid->x[p] = x/Rho;
    solid->y[p] = y/Rho;
    solid->z[p] = z/Rho*(Dim-2);

  }

  return 0;
}

PetscErrorCode SoftForceOnOneEllipsoid(System_Struct *system, Solid_Struct *solid)
{
  PetscScalar E, dEdx, dEdy, dEdz;
  const PetscInt Dim = system->Dim;
  const PetscScalar A  = solid->forces.SF.A;
  const PetscScalar B  = solid->forces.SF.B;
  const PetscScalar C  = solid->forces.SF.C;
  const PetscScalar sf = solid->forces.SF.sf;

  for (PetscInt j=0; j<solid->forces.SF.N; j++)
  {
    const PetscInt p = solid->forces.SF.p[j];

    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    E = -1.0 + x*x/(A*A) + y*y/(B*B) + (Dim-2)*z*z/(C*C); 

    
    dEdx = 2.0*x/(A*A);
    dEdy = 2.0*y/(B*B);
    dEdz = 2.0*z/(C*C);

    solid->forces.Fx[p] += -sf*E*dEdx;
    solid->forces.Fy[p] += -sf*E*dEdy;
    solid->forces.Fz[p] += -sf*E*dEdz;

  }
  return 0;
}

PetscErrorCode SoftForceOnTwoEllipsoids(System_Struct *system, Solid_Struct *solid)
{
  PetscScalar E, dEdx, dEdy, dEdz;
  const PetscInt Dim = system->Dim;

  const PetscScalar A  = solid->forces.SF.A;
  const PetscScalar Bd = solid->forces.SF.Bd;
  const PetscScalar Bv = solid->forces.SF.Bv;
  const PetscScalar C  = solid->forces.SF.C;
  const PetscScalar sf = solid->forces.SF.sf;
  PetscScalar B;

  for (PetscInt j=0; j<solid->forces.SF.N; j++)
  {
    const PetscInt p = solid->forces.SF.p[j];

    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    if (y>0) B = Bd;
    else     B = Bv;
    
    E = -1.0 + x*x/(A*A) + y*y/(B*B) + (Dim-2)*z*z/(C*C); 

    
    dEdx = 2.0*x/(A*A);
    dEdy = 2.0*y/(B*B);
    dEdz = 2.0*z/(C*C);

    solid->forces.Fx[p] += -sf*E*dEdx;
    solid->forces.Fy[p] += -sf*E*dEdy;
    solid->forces.Fz[p] += -sf*E*dEdz;

  }
  return 0;
}


PetscErrorCode SoftForceNull(System_Struct *system, Solid_Struct *solid) { return 0;}

PetscErrorCode AddSoftForce(System_Struct *system, Solid_Struct *solid)
{
  if (solid->forces.SF.Type == FORCE_SF_ELLIPSOID)
    SoftForceOnOneEllipsoid(system, solid);
  else if (solid->forces.SF.Type == FORCE_SF_TWO_ELLIPSOIDS)
    SoftForceOnTwoEllipsoids(system, solid);
  return 0;
}

PetscErrorCode SoftForceInitialization(System_Struct *system, Solid_Struct *solid)
{   
  if (solid->forces.SF.Type == FORCE_SF_ELLIPSOID)
    SoftForceCreateOneEllipsoid(system, solid);
  else if (solid->forces.SF.Type == FORCE_SF_TWO_ELLIPSOIDS)
    SoftForceCreateTwoEllipsoids(system, solid);
  return 0;
}



#endif