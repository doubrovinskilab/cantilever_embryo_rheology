#ifndef ElementalMatrixBoundaryRobin_C
#define ElementalMatrixBoundaryRobin_C

PetscErrorCode Calc2DP2_Urlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  
  PetscInt **e_bound, PetscScalar **Urlp)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar Dx = x1 - x0, Dy = y1-y0;
  const PetscScalar L = PetscSqrtScalar( Dx*Dx + Dy*Dy );

  Urlp[0][0] =  2.0*L/15.0;   Urlp[0][1] =     -L/30.0;  Urlp[0][2] =      L/15.0;
  Urlp[1][0] =     -L/30.0;   Urlp[1][1] =  2.0*L/15.0;  Urlp[1][2] =      L/15.0;
  Urlp[2][0] =      L/15.0;   Urlp[2][1] =      L/15.0;  Urlp[2][2] =  8.0*L/15.0;

  return 0;
}

PetscErrorCode Calc2DP1b_Urlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  
  PetscInt **e_bound, PetscScalar **Urlp)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  
  const PetscScalar Dx = x1 - x0, Dy = y1-y0;
  const PetscScalar L = PetscSqrtScalar( Dx*Dx + Dy*Dy );

  Urlp[0][0] =  L/3.0;   Urlp[0][1] =  L/6.0;
  Urlp[1][0] =  L/6.0;   Urlp[1][1] =  L/3.0;
  
  return 0;
}

PetscErrorCode Calc3DP2_Urlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  
  PetscInt **e_bound, PetscScalar **Urlp)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscInt n2 = e_bound[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];

  const PetscScalar Dx10 = x1-x0, Dy10 = y1-y0, Dz10 = z1-z0;
  const PetscScalar Dx20 = x2-x0, Dy20 = y2-y0, Dz20 = z2-z0;

  const PetscScalar Area = 0.5*PetscSqrtScalar( PetscPowScalar((Dx10*Dy20 - Dx20*Dy10),2) 
    + PetscPowScalar((-Dx10*Dz20 + Dx20*Dz10),2) + PetscPowScalar((Dy10*Dz20 - Dy20*Dz10),2) );

  Urlp[0][0] =     Area/30.0;
  Urlp[0][1] =   -Area/180.0;
  Urlp[0][2] =   -Area/180.0;
  Urlp[0][3] =           0.0;
  Urlp[0][4] =    -Area/45.0;
  Urlp[0][5] =           0.0;
  Urlp[1][0] =   -Area/180.0;
  Urlp[1][1] =     Area/30.0;
  Urlp[1][2] =   -Area/180.0;
  Urlp[1][3] =           0.0;
  Urlp[1][4] =           0.0;
  Urlp[1][5] =    -Area/45.0;
  Urlp[2][0] =   -Area/180.0;
  Urlp[2][1] =   -Area/180.0;
  Urlp[2][2] =     Area/30.0;
  Urlp[2][3] =    -Area/45.0;
  Urlp[2][4] =           0.0;
  Urlp[2][5] =           0.0;
  Urlp[3][0] =           0.0;
  Urlp[3][1] =           0.0;
  Urlp[3][2] =    -Area/45.0;
  Urlp[3][3] = 8.0*Area/45.0;
  Urlp[3][4] = 4.0*Area/45.0;
  Urlp[3][5] = 4.0*Area/45.0;
  Urlp[4][0] =    -Area/45.0;
  Urlp[4][1] =           0.0;
  Urlp[4][2] =           0.0;
  Urlp[4][3] = 4.0*Area/45.0;
  Urlp[4][4] = 8.0*Area/45.0;
  Urlp[4][5] = 4.0*Area/45.0;
  Urlp[5][0] =           0.0;
  Urlp[5][1] =    -Area/45.0;
  Urlp[5][2] =           0.0;
  Urlp[5][3] = 4.0*Area/45.0;
  Urlp[5][4] = 4.0*Area/45.0;
  Urlp[5][5] = 8.0*Area/45.0;

  return 0;
}

PetscErrorCode Calc3DP1b_Urlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  
  PetscInt **e_bound, PetscScalar **Urlp)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscInt n2 = e_bound[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];

  const PetscScalar Dx10 = x1-x0, Dy10 = y1-y0, Dz10 = z1-z0;
  const PetscScalar Dx20 = x2-x0, Dy20 = y2-y0, Dz20 = z2-z0;

  const PetscScalar Area = 0.5*PetscSqrtScalar( PetscPowScalar((Dx10*Dy20 - Dx20*Dy10),2) 
    + PetscPowScalar((-Dx10*Dz20 + Dx20*Dz10),2) + PetscPowScalar((Dy10*Dz20 - Dy20*Dz10),2) );

  Urlp[0][0] =  Area/6.;
  Urlp[0][1] =  Area/12.;
  Urlp[0][2] =  Area/12.;

  Urlp[1][0] =  Area/12.;
  Urlp[1][1] =  Area/6.;
  Urlp[1][2] =  Area/12.;

  Urlp[2][0] =  Area/12.;
  Urlp[2][1] =  Area/12.;
  Urlp[2][2] =  Area/6.;

  return 0;
}


#endif