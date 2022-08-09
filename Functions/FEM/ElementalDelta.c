#ifndef ElementalDelta_C
#define ElementalDelta_C

PetscErrorCode Calc2D_Delta( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, PetscScalar *Delta)
{
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];

  const PetscScalar x0 = _x[n0], y0 = _y[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2];

  *Delta =  x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1; // Delta

  return 0;
}

PetscErrorCode Calc3D_Delta( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z,
  PetscInt **_inter, PetscScalar *Delta)
{  
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];
  const PetscInt n3 = _inter[en][3];

  const PetscScalar x0 = _x[n0], y0 = _y[n0], z0 = _z[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1], z1 = _z[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2], z2 = _z[n2];
  const PetscScalar x3 = _x[n3], y3 = _y[n3], z3 = _z[n3];

  *Delta =  -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2  
    + x1*y0*z2 - x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3 
    + x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0 + x3*y1*z2 
    + x3*y2*z0 - x3*y2*z1;

    return 0;
}

#endif