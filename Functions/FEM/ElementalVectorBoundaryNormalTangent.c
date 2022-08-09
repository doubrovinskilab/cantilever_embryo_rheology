#ifndef ElementalVectorBoundaryNormalTangent_C
#define ElementalVectorBoundaryNormalTangent_C


PetscErrorCode Calc2D_NormTang(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_bound, PetscScalar *N, 
  PetscScalar *T1, PetscScalar *T2, const PetscScalar *DC)
{

  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];

  const PetscScalar x0 = n_x[n0],     y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1],     y1 = n_y[n1];
  const PetscScalar xc = (x0+x1)/2.,  yc = (y0+y1)/2.;
  
  PetscScalar Dx = x1-x0,  Dy = y1-y0;
  PetscScalar L = PetscSqrtScalar(Dx*Dx + Dy*Dy);
  T1[0] = Dx / L;
  T1[1] = Dy / L;
  T1[2] = 0.0;

  PetscScalar Nx, Ny;
  if (Dx != 0.0)  Nx = -Dy/Dx, Ny = 1.0; 
  else            Nx = -1.0, Ny = Dx/Dy; 

  const PetscScalar Cx = xc - DC[0], Cy = yc - DC[1]; 
  const PetscScalar DOT = Cx*Nx + Cy*Ny;              

  if (DOT < 0) 
    Nx = -1*Nx, Ny = -1*Ny;
  L = PetscSqrtScalar(Nx*Nx + Ny*Ny);
  N[0] = Nx/L;
  N[1] = Ny/L;
  N[2] = 0.0;

  T2[0] = 0;
  T2[1] = 0;
  T2[2] = 0;

  return 0;
}


PetscErrorCode Calc3D_NormTang(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,PetscInt **e_bound, PetscScalar *N, 
  PetscScalar *T1, PetscScalar *T2, const PetscScalar *DC)
{

  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscInt n2 = e_bound[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar xc = (x0 + x1 + x2)/3.0;
  const PetscScalar yc = (y0 + y1 + y2)/3.0;
  const PetscScalar zc = (z0 + z1 + z2)/3.0;
  
  const PetscScalar Dx10 = x1-x0, Dy10 = y1-y0, Dz10 = z1-z0;
  const PetscScalar Dx20 = x2-x0, Dy20 = y2-y0, Dz20 = z2-z0;
  CalcCrossProduct( Dx10, Dy10, Dz10, Dx20, Dy20, Dz20, &N[0], &N[1], &N[2]);

  const PetscScalar Cx = xc - DC[0], Cy = yc - DC[1], Cz = zc - DC[2]; 
  const PetscScalar DOT = Cx*N[0] + Cy*N[1] + Cz*N[2];              

  if (DOT < 0) 
    N[0] = -1.0*N[0], N[1] = -1.0*N[1], N[2] = -1.0*N[2];
  PetscScalar L = PetscSqrtScalar(N[0]*N[0] + N[1]*N[1] + N[2]*N[2]);
  N[0] = N[0]/L;
  N[1] = N[1]/L;
  N[2] = N[2]/L;

  
  L = PetscSqrtScalar(Dx20*Dx20 + Dy20*Dy20 + Dz20*Dz20);
  T2[0] = Dx20/L;
  T2[1] = Dy20/L;
  T2[2] = Dz20/L;

  
  CalcCrossProduct( T2[0], T2[1], T2[0], N[0], N[1], N[2], &T1[0], &T1[1], &T1[2]);
  L = PetscSqrtScalar(T1[0]*T1[0] + T1[1]*T1[1] + T1[2]*T1[2]);
  T1[0] = T1[0]/L;
  T1[1] = T1[1]/L;
  T1[2] = T1[2]/L;

  return 0;
}
#endif