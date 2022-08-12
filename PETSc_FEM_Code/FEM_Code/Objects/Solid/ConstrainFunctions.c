#ifndef Solid_Constrain_Functions_C
#define Solid_Constrain_Functions_C

PetscErrorCode ConstrainCorrectPointsOnOneEllipsoid(System_Struct *system, Solid_Struct *solid, PetscInt *ps_mpi, 
  PetscScalar *xs_mpi, PetscScalar *ys_mpi, PetscScalar *zs_mpi);
PetscErrorCode ConstrainCorrectPointsOnOneEllipsoidSerial(System_Struct *system, Solid_Struct *solid);
PetscErrorCode ConstrainCorrectPointsOnTwoEllipsoids(System_Struct *system, Solid_Struct *solid, PetscInt *ps_mpi, 
  PetscScalar *xs_mpi, PetscScalar *ys_mpi, PetscScalar *zs_mpi);
PetscErrorCode ConstrainCorrectPointsOnTwoEllipsoidsSerial(System_Struct *system, Solid_Struct *solid);
PetscErrorCode ConstrainCorrectPointsOnAxis(System_Struct *system, Solid_Struct *solid, PetscInt *ps_mpi, 
  PetscScalar *xs_mpi, PetscScalar *ys_mpi, PetscScalar *zs_mpi);
PetscErrorCode ConstrainCorrectPointsOnAxisSerial(System_Struct *system, Solid_Struct *solid);

PetscErrorCode ConstrainCreateOneEllipsoid(System_Struct *system, Solid_Struct *solid)
{
  const PetscInt Dim = system->Dim;

  const PetscScalar ds = solid->constrain.ds;
  const PetscScalar As = solid->constrain.A - ds;
  const PetscScalar Al = solid->constrain.A + ds;
  const PetscScalar Bs = solid->constrain.B - ds;
  const PetscScalar Bl = solid->constrain.B + ds;
  const PetscScalar Cs = solid->constrain.C - ds;
  const PetscScalar Cl = solid->constrain.C + ds;

  PetscScalar Es, El = 0; 
  
  solid->constrain.N = 0;
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

    if (Es>=1.0 && El<=1.0) solid->constrain.N++;
  }

  PetscMalloc1(solid->constrain.N, &solid->constrain.p);

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
      solid->constrain.p[j] = p;
      j++;}
  }
  ConstrainCorrectPointsOnOneEllipsoidSerial(system, solid); 
  return 0;
}

PetscErrorCode ConstrainCorrectPointsOnOneEllipsoid(System_Struct *system, Solid_Struct *solid, PetscInt *ps_mpi, 
  PetscScalar *xs_mpi, PetscScalar *ys_mpi, PetscScalar *zs_mpi)
{
  PetscScalar Rho;
  PetscBool found_p;

  const PetscInt Dim = system->Dim;

  const PetscScalar A  = solid->constrain.A;
  const PetscScalar B  = solid->constrain.B;
  const PetscScalar C  = solid->constrain.C;

  for (PetscInt j=0; j<solid->constrain.N; j++)
  {
    const PetscInt p = solid->constrain.p[j];

    if (solid->fluid_elem[p] == -1)
      continue; 

    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    Rho = PetscSqrtScalar( (x*x)/(A*A) + (y*y)/(B*B) + (z*z)/(C*C) );

    solid->x[p] = x/Rho;
    solid->y[p] = y/Rho;
    solid->z[p] = z/Rho*(Dim-2);

    
    PetscInt i = 0;
    found_p = PETSC_FALSE;
    while(PETSC_TRUE) 
    {
      if( ps_mpi[i] == p) {
        found_p = PETSC_TRUE;
        break; }
      i++; 
    }

    if (found_p)
    { xs_mpi[i] = solid->x[p], ys_mpi[i] = solid->y[p], zs_mpi[i] = solid->z[p];}
  }
  return 0;
}


PetscErrorCode ConstrainCorrectPointsOnOneEllipsoidSerial(System_Struct *system, Solid_Struct *solid)
{
  PetscScalar Rho;

  const PetscInt Dim = system->Dim;

  const PetscScalar A  = solid->constrain.A;
  const PetscScalar B  = solid->constrain.B;
  const PetscScalar C  = solid->constrain.C;

  for (PetscInt j=0; j<solid->constrain.N; j++)
  {
    const PetscInt p = solid->constrain.p[j];

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

PetscErrorCode ConstrainCreateTwoEllipsoids(System_Struct *system, Solid_Struct *solid)
{
  const PetscInt Dim = system->Dim;

  const PetscScalar ds  = solid->constrain.ds;
  const PetscScalar As  = solid->constrain.A  - ds;
  const PetscScalar Al  = solid->constrain.A  + ds;
  const PetscScalar Bds = solid->constrain.Bd - ds;
  const PetscScalar Bdl = solid->constrain.Bd + ds;
  const PetscScalar Bvs = solid->constrain.Bv - ds;
  const PetscScalar Bvl = solid->constrain.Bv + ds;
  const PetscScalar Cs  = solid->constrain.C  - ds;
  const PetscScalar Cl  = solid->constrain.C  + ds;
  PetscScalar Bs, Bl;
  PetscScalar Es, El = 0; 

  solid->constrain.N = 0;
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

    if (Es>=1.0 && El<=1.0) solid->constrain.N++;
  }

  PetscMalloc1(solid->constrain.N, &solid->constrain.p);

  // Saving the points on the constrain
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
      solid->constrain.p[j] = p;
      j++;}
  }

  ConstrainCorrectPointsOnTwoEllipsoidsSerial(system, solid); 

  return 0;
}

PetscErrorCode ConstrainCorrectPointsOnTwoEllipsoids(System_Struct *system, Solid_Struct *solid, PetscInt *ps_mpi, 
  PetscScalar *xs_mpi, PetscScalar *ys_mpi, PetscScalar *zs_mpi)
{
  PetscScalar Rho;
  PetscBool found_p;

  const PetscInt Dim = system->Dim;

  const PetscScalar A  = solid->constrain.A;
  const PetscScalar Bd = solid->constrain.Bd;
  const PetscScalar Bv = solid->constrain.Bv;
  const PetscScalar C  = solid->constrain.C;
  PetscScalar B;

  for (PetscInt j=0; j<solid->constrain.N; j++)
  {
    const PetscInt p = solid->constrain.p[j];

    if (solid->fluid_elem[p] == -1)
      continue; 

    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    if (y>0) B = Bd;
    else     B = Bv;

    Rho = PetscSqrtScalar( (x*x)/(A*A) + (y*y)/(B*B) + (z*z)/(C*C) );

    solid->x[p] = x/Rho;
    solid->y[p] = y/Rho;
    solid->z[p] = z/Rho*(Dim-2);

    PetscInt i = 0;
    found_p = PETSC_FALSE;
    while(PETSC_TRUE) 
    {
      if( ps_mpi[i] == p) {
        found_p = PETSC_TRUE;
        break; }
      i++; 
    }

    if (found_p)
    { xs_mpi[i] = solid->x[p], ys_mpi[i] = solid->y[p], zs_mpi[i] = solid->z[p];}
  }

  return 0;
}


PetscErrorCode ConstrainCorrectPointsOnTwoEllipsoidsSerial(System_Struct *system, Solid_Struct *solid)
{
  PetscScalar Rho;
  const PetscInt Dim = system->Dim;

  const PetscScalar A  = solid->constrain.A;
  const PetscScalar Bd = solid->constrain.Bd;
  const PetscScalar Bv = solid->constrain.Bv;
  const PetscScalar C  = solid->constrain.C;
  PetscScalar B;

  for (PetscInt j=0; j<solid->constrain.N; j++)
  {
    const PetscInt p = solid->constrain.p[j];

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


PetscErrorCode ConstrainCreateAxis(System_Struct *system, Solid_Struct *solid)
{
  const PetscScalar ds    = solid->constrain.ds;
  const PetscScalar Vmin  = solid->constrain.axis_v - ds;
  const PetscScalar Vmax  = solid->constrain.axis_v + ds;
  const PetscInt    axis  = solid->constrain.axis_i;

  if (axis!=0 && axis==1 && axis==2)
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in SolidCreateAxisConstrain");  

  solid->constrain.N = 0;
  for (PetscInt p=0; p<solid->N_total_pts; p++)
  {
    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];
    if (axis==0 && x>=Vmin && x<=Vmax)
        solid->constrain.N++;
    else if (axis==1 && y>=Vmin && y<=Vmax)
        solid->constrain.N++;
    else if (axis==2 && z>=Vmin && z<=Vmax)
        solid->constrain.N++;
  }
  PetscMalloc1(solid->constrain.N, &solid->constrain.p);

  PetscInt j = 0;
  for (PetscInt p=0; p<solid->N_total_pts; p++)
  {
    const PetscScalar x = solid->x[p];
    const PetscScalar y = solid->y[p];
    const PetscScalar z = solid->z[p];

    if (axis==0 && x>=Vmin && x<=Vmax){
      solid->constrain.p[j] = p;    j++; }
    else if (axis==1 && y>=Vmin && y<=Vmax){
      solid->constrain.p[j] = p;    j++; }
    else if (axis==2 && z>=Vmin && z<=Vmax){
      solid->constrain.p[j] = p;    j++; }
  }
  PetscPrintf(PETSC_COMM_WORLD,"    The number of constrained points is %d\n",  solid->constrain.N); 
  ConstrainCorrectPointsOnAxisSerial(system, solid);

  return 0;
}

PetscErrorCode ConstrainCorrectPointsOnAxis(System_Struct *system, Solid_Struct *solid, PetscInt *ps_mpi, 
  PetscScalar *xs_mpi, PetscScalar *ys_mpi, PetscScalar *zs_mpi)
{
  PetscBool found_p;

  const PetscScalar V = solid->constrain.axis_v;
  const PetscInt axis = solid->constrain.axis_i;
  const PetscInt N    = solid->constrain.N;

  for (PetscInt j=0; j<N; j++)
  {
    const PetscInt p = solid->constrain.p[j];

    if (solid->fluid_elem[p] == -1)
      continue;

    if (axis==0)       solid->x[p] = V;
    else if  (axis==1) solid->y[p] = V;
    else if  (axis==2) solid->z[p] = V;

    // Adjust Data to be sent
    PetscInt i = 0;
    found_p = PETSC_FALSE;
    while(PETSC_TRUE) 
    {
      if( ps_mpi[i] == p) {
        found_p = PETSC_TRUE;
        break; }
      i++; 
    }

    if (found_p)
    { xs_mpi[i] = solid->x[p], ys_mpi[i] = solid->y[p], zs_mpi[i] = solid->z[p];}
  }
  return 0;
}


PetscErrorCode ConstrainCorrectPointsOnAxisSerial(System_Struct *system, Solid_Struct *solid)
{
  const PetscScalar V = solid->constrain.axis_v;
  const PetscInt axis = solid->constrain.axis_i;
  const PetscInt N    = solid->constrain.N;

  for (PetscInt j=0; j<N; j++)
  {
    const PetscInt p = solid->constrain.p[j];

    if (axis==0)       solid->x[p] = V;
    else if  (axis==1) solid->y[p] = V;
    else if  (axis==2) solid->z[p] = V;
  }
  return 0;
}

PetscErrorCode ConstrainCreate(System_Struct *system, Solid_Struct *solid)
{
  if (solid->constrain.Type == CONSTRAIN_ELLIPSOID)
    ConstrainCreateOneEllipsoid(system, solid);
  else if (solid->constrain.Type == CONSTRAIN_TWO_ELLIPSOIDS)
    ConstrainCreateTwoEllipsoids(system, solid);
  else if (solid->constrain.Type == CONSTRAIN_AXIS)
    ConstrainCreateAxis(system, solid);
  else
    return 0;

  return 0;
}

PetscErrorCode ConstrainCorrectPoints(System_Struct *system, Solid_Struct *solid, PetscInt *ps_mpi, 
  PetscScalar *xs_mpi, PetscScalar *ys_mpi, PetscScalar *zs_mpi)
{
  if (solid->constrain.Type == CONSTRAIN_ELLIPSOID)
    ConstrainCorrectPointsOnOneEllipsoid(system, solid, ps_mpi, xs_mpi, ys_mpi, zs_mpi);
  else if (solid->constrain.Type == CONSTRAIN_ELLIPSOID)
    ConstrainCorrectPointsOnTwoEllipsoids(system, solid, ps_mpi, xs_mpi, ys_mpi, zs_mpi);
  else if (solid->constrain.Type == CONSTRAIN_AXIS)
    ConstrainCorrectPointsOnAxis(system, solid, ps_mpi, xs_mpi, ys_mpi, zs_mpi);
  else
    return 0;
  return 0;
}

PetscErrorCode ConstrainCorrectPointsSerial(System_Struct *system, Solid_Struct *solid)
{
  if (solid->constrain.Type == CONSTRAIN_ELLIPSOID)
    ConstrainCorrectPointsOnOneEllipsoidSerial(system, solid);
  else if (solid->constrain.Type == CONSTRAIN_ELLIPSOID)
    ConstrainCorrectPointsOnTwoEllipsoidsSerial(system, solid);
  else if (solid->constrain.Type == CONSTRAIN_AXIS)
    ConstrainCorrectPointsOnAxisSerial(system, solid);
  else
    return 0;
  return 0;
}

PetscErrorCode ConstrainFree(Solid_Struct *solid)
{
  if (solid->constrain.Type == CONSTRAIN_NULL)
    PetscFree(solid->constrain.p);
  return 0;
}

#endif 
