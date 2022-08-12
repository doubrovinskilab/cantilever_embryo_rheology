#ifndef ElementalVectorBoundaryNeumann_C
#define ElementalVectorBoundaryNeumann_C

PetscErrorCode Calc2DP1_Unlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, PetscScalar *Unlp)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  
  const PetscScalar Dx = x1-x0, Dy = y1-y0;
  const PetscScalar L = PetscSqrtScalar(Dx*Dx + Dy*Dy);

  Unlp[0] = 0.5*L;
  Unlp[1] = 0.5*L;
  return 0;
}


PetscErrorCode Calc2DP2_Unlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, PetscScalar *Unlp)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  
  const PetscScalar Dx = x1-x0, Dy = y1-y0;
  const PetscScalar L = PetscSqrtScalar(Dx*Dx + Dy*Dy);

  Unlp[0] =     L/6.0;
  Unlp[1] =     L/6.0;
  Unlp[2] = 4.0*L/6.0;
  return 0;
}

PetscErrorCode Calc2DP1b_Unlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, PetscScalar *Unlp)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  
  const PetscScalar Dx = x1-x0, Dy = y1-y0;
  const PetscScalar L = PetscSqrtScalar(Dx*Dx + Dy*Dy);

  Unlp[0] =     L/2.0;
  Unlp[1] =     L/2.0;
  return 0;
}


PetscErrorCode Calc3DP1_Unlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, PetscScalar *Unlp)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscInt n2 = e_bound[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  
  const PetscScalar Dx10 = x1-x0, Dy10 = y1-y0, Dz10 =z1-z0;
  const PetscScalar Dx20 = x2-x0, Dy20 = y2-y0, Dz20 =z2-z0;

  const PetscScalar Area = 0.5*PetscSqrtScalar( PetscPowScalar((Dx10*Dy20 - Dx20*Dy10),2) + PetscPowScalar((-Dx10*Dz20 + Dx20*Dz10),2) + PetscPowScalar((Dy10*Dz20 - Dy20*Dz10),2) );

  Unlp[0] = Area/3.0;
  Unlp[1] = Area/3.0;
  Unlp[2] = Area/3.0;

  return 0;
}

PetscErrorCode Calc3DP2_Unlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, PetscScalar *Unlp)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscInt n2 = e_bound[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  
  const PetscScalar Dx10 = x1-x0, Dy10 = y1-y0, Dz10 =z1-z0;
  const PetscScalar Dx20 = x2-x0, Dy20 = y2-y0, Dz20 =z2-z0;

  const PetscScalar Area = 0.5*PetscSqrtScalar( PetscPowScalar((Dx10*Dy20 - Dx20*Dy10),2) + PetscPowScalar((-Dx10*Dz20 + Dx20*Dz10),2) + PetscPowScalar((Dy10*Dz20 - Dy20*Dz10),2) );

  Unlp[0] = 0.;
  Unlp[1] = 0.;
  Unlp[2] = 0.;
  Unlp[3] = Area/3.;
  Unlp[4] = Area/3.;
  Unlp[5] = Area/3.;

  return 0;
}

PetscErrorCode Calc3DP1b_Unlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, PetscScalar *Unlp)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscInt n2 = e_bound[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  
  const PetscScalar Dx10 = x1-x0, Dy10 = y1-y0, Dz10 =z1-z0;
  const PetscScalar Dx20 = x2-x0, Dy20 = y2-y0, Dz20 =z2-z0;

  const PetscScalar Area = 0.5*PetscSqrtScalar( PetscPowScalar((Dx10*Dy20 - Dx20*Dy10),2) + PetscPowScalar((-Dx10*Dz20 + Dx20*Dz10),2) + PetscPowScalar((Dy10*Dz20 - Dy20*Dz10),2) );

  Unlp[0] = Area/3.;
  Unlp[1] = Area/3.;
  Unlp[2] = Area/3.;

  return 0;
}
#endif