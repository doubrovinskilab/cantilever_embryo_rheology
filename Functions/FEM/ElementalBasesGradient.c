#ifndef ElementaBasesGradient_C
#define ElementaBasesGradient_C

PetscErrorCode Calc2DP1_GradPhi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z,
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp,
  PetscScalar *grad_phi_x, PetscScalar *grad_phi_y, PetscScalar *grad_phi_z)
{
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];

  const PetscScalar x0 = _x[n0], y0 = _y[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2];

  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;

  grad_phi_x[0] =  (y1 - y2)/Delta;
  grad_phi_x[1] = -(y0 - y2)/Delta;
  grad_phi_x[2] =  (y0 - y1)/Delta;
  
  grad_phi_y[0] = -(x1 - x2)/Delta;
  grad_phi_y[1] =  (x0 - x2)/Delta;
  grad_phi_y[2] = -(x0 - x1)/Delta;

  grad_phi_z[0] = 0;
  grad_phi_z[1] = 0;
  grad_phi_z[2] = 0;

  return 0;
}


PetscErrorCode Calc2DP1b_GradPhi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z,
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp,
  PetscScalar *grad_phi_x, PetscScalar *grad_phi_y, PetscScalar *grad_phi_z)
{
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];

  const PetscScalar x0 = _x[n0], y0 = _y[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2];

  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;
  const PetscScalar Delta3 = Delta*Delta*Delta;

  grad_phi_x[0] =  (y1 - y2)/Delta + 9*(-y0 + y1)*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(-y0 + y2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(y1 - y2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))/Delta3;
  grad_phi_x[1] =  (-y0 + y2)/Delta + 9*(-y0 + y1)*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(-y0 + y2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(y1 - y2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))/Delta3;
  grad_phi_x[2] = -(-y0 + y1)/Delta + 9*(-y0 + y1)*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(-y0 + y2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(y1 - y2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))/Delta3;
  grad_phi_x[3] = -27*(-y0 + y1)*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 - 27*(-y0 + y2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 - 27*(y1 - y2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))/Delta3;
  
  grad_phi_y[0] =  (-x1 + x2)/Delta + 9*(x0 - x1)*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(x0 - x2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(-x1 + x2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))/Delta3;
  grad_phi_y[1] =   (x0 - x2)/Delta + 9*(x0 - x1)*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(x0 - x2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(-x1 + x2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))/Delta3;
  grad_phi_y[2] =  -(x0 - x1)/Delta + 9*(x0 - x1)*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(x0 - x2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 + 9*(-x1 + x2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))/Delta3;
  grad_phi_y[3] = -27*(x0 - x1)*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 - 27*(x0 - x2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(y1 - y2) + x1*y2 - x2*y1 + yp*(-x1 + x2))/Delta3 - 27*(-x1 + x2)*(xp*(-y0 + y1) - x0*y1 + x1*y0 + yp*(x0 - x1))*(xp*(-y0 + y2) - x0*y2 + x2*y0 + yp*(x0 - x2))/Delta3;

  grad_phi_z[0] = 0;
  grad_phi_z[1] = 0;
  grad_phi_z[2] = 0;
  grad_phi_z[3] = 0;

  return 0;
}


PetscErrorCode Calc2DP2_GradPhi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z,
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp,
  PetscScalar *grad_phi_x, PetscScalar *grad_phi_y, PetscScalar *grad_phi_z)
{   
  const PetscInt n0  = _inter[en][0];
  const PetscInt n1  = _inter[en][1];
  const PetscInt n2  = _inter[en][2];

  const PetscScalar x0 = _x[n0], y0 = _y[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2];

  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;
  const PetscScalar Delta2 = Delta*Delta;

  grad_phi_x[0] =  (y1 - y2)*(-Delta + 4*xp*y1 - 4*xp*y2 - 4*x1*yp + 4*x1*y2 + 4*x2*yp - 4*x2*y1)/Delta2;
  grad_phi_x[1] = -(y0 - y2)*(-Delta - 4*xp*y0 + 4*xp*y2 + 4*x0*yp - 4*x0*y2 - 4*x2*yp + 4*x2*y0)/Delta2;
  grad_phi_x[2] = -(y0 - y1)*( Delta - 4*xp*y0 + 4*xp*y1 + 4*x0*yp - 4*x0*y1 - 4*x1*yp + 4*x1*y0)/Delta2,
  grad_phi_x[3] =  4*(-2*xp*y0*y1 + 2*xp*y0*y2 + 2*xp*y1*y2 - 2*xp*y2*y2 + x0*yp*y1 - x0*yp*y2 - x0*y1*y2 + x0*y2*y2 
                        + x1*yp*y0 - x1*yp*y2 - x1*y0*y2 + x1*y2*y2 - x2*yp*y0 - x2*yp*y1 + 2*x2*yp*y2 + 2*x2*y0*y1 - x2*y0*y2 - x2*y1*y2)/Delta2;
  grad_phi_x[4] =  4*(-2*xp*y0*y0 + 2*xp*y0*y1 + 2*xp*y0*y2 - 2*xp*y1*y2 + 2*x0*yp*y0 - x0*yp*y1 - x0*yp*y2 - x0*y0*y1 
                        - x0*y0*y2 + 2*x0*y1*y2 - x1*yp*y0 + x1*yp*y2 + x1*y0*y0 - x1*y0*y2 - x2*yp*y0 + x2*yp*y1 + x2*y0*y0 - x2*y0*y1)/Delta2;
  grad_phi_x[5] = -4*(-2*xp*y0*y1 + 2*xp*y0*y2 + 2*xp*y1*y1 - 2*xp*y1*y2 + x0*yp*y1 - x0*yp*y2 - x0*y1*y1 + x0*y1*y2 
                        + x1*yp*y0 - 2*x1*yp*y1 + x1*yp*y2 + x1*y0*y1 - 2*x1*y0*y2 + x1*y1*y2 - x2*yp*y0 + x2*yp*y1 + x2*y0*y1 - x2*y1*y1)/Delta2;

  grad_phi_y[0] =  -(x1 - x2)*(-Delta + 4*xp*y1 - 4*xp*y2 - 4*x1*yp + 4*x1*y2 + 4*x2*yp - 4*x2*y1)/Delta2;
  grad_phi_y[1] =   (x0 - x2)*(-Delta - 4*xp*y0 + 4*xp*y2 + 4*x0*yp - 4*x0*y2 - 4*x2*yp + 4*x2*y0)/Delta2;
  grad_phi_y[2] =   (x0 - x1)*( Delta - 4*xp*y0 + 4*xp*y1 + 4*x0*yp - 4*x0*y1 - 4*x1*yp + 4*x1*y0)/Delta2;
  grad_phi_y[3] =   4*(xp*x0*y1 - xp*x0*y2 + xp*x1*y0 - xp*x1*y2 - xp*x2*y0 - xp*x2*y1 + 2*xp*x2*y2 - 2*x0*x1*yp + 2*x0*x1*y2 
                         + 2*x0*x2*yp - x0*x2*y1 - x0*x2*y2 + 2*x1*x2*yp - x1*x2*y0 - x1*x2*y2 - 2*x2*x2*yp + x2*x2*y0 + x2*x2*y1)/Delta2;
  grad_phi_y[4] =  -4*(-2*xp*x0*y0 + xp*x0*y1 + xp*x0*y2 + xp*x1*y0 - xp*x1*y2 + xp*x2*y0 - xp*x2*y1 + 2*x0*x0*yp - x0*x0*y1 
                         - x0*x0*y2 - 2*x0*x1*yp + x0*x1*y0 + x0*x1*y2 - 2*x0*x2*yp + x0*x2*y0 + x0*x2*y1 + 2*x1*x2*yp - 2*x1*x2*y0)/Delta2;
  grad_phi_y[5] =  -4*(xp*x0*y1 - xp*x0*y2 + xp*x1*y0 - 2*xp*x1*y1 + xp*x1*y2 - xp*x2*y0 + xp*x2*y1 - 2*x0*x1*yp + x0*x1*y1 
                         + x0*x1*y2 + 2*x0*x2*yp - 2*x0*x2*y1 + 2*x1*x1*yp - x1*x1*y0 - x1*x1*y2 - 2*x1*x2*yp + x1*x2*y0 + x1*x2*y1)/Delta2;

  grad_phi_z[0] = 0;
  grad_phi_z[1] = 0;
  grad_phi_z[2] = 0;
  grad_phi_z[3] = 0;
  grad_phi_z[4] = 0;
  grad_phi_z[5] = 0;

  return 0;
}

PetscErrorCode Calc3DP1_GradPhi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp,
  PetscScalar *grad_phi_x, PetscScalar *grad_phi_y, PetscScalar *grad_phi_z)
{  
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];
  const PetscInt n3 = _inter[en][3];

  const PetscScalar x0 = _x[n0], y0 = _y[n0], z0 = _z[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1], z1 = _z[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2], z2 = _z[n2];
  const PetscScalar x3 = _x[n3], y3 = _y[n3], z3 = _z[n3];

  const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2  
    + x1*y0*z2 - x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3 
    + x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0 + x3*y1*z2 
    + x3*y2*z0 - x3*y2*z1;
  
  grad_phi_x[0] = -(y1*z2 - y1*z3 - y2*z1 + y2*z3 + y3*z1 - y3*z2)/Delta;
  grad_phi_x[1] =  (y0*z2 - y0*z3 - y2*z0 + y2*z3 + y3*z0 - y3*z2)/Delta;
  grad_phi_x[2] = -(y0*z1 - y0*z3 - y1*z0 + y1*z3 + y3*z0 - y3*z1)/Delta;
  grad_phi_x[3] =  (y0*z1 - y0*z2 - y1*z0 + y1*z2 + y2*z0 - y2*z1)/Delta;

  grad_phi_y[0] =  (x1*z2 - x1*z3 - x2*z1 + x2*z3 + x3*z1 - x3*z2)/Delta;
  grad_phi_y[1] = -(x0*z2 - x0*z3 - x2*z0 + x2*z3 + x3*z0 - x3*z2)/Delta;
  grad_phi_y[2] =  (x0*z1 - x0*z3 - x1*z0 + x1*z3 + x3*z0 - x3*z1)/Delta;
  grad_phi_y[3] = -(x0*z1 - x0*z2 - x1*z0 + x1*z2 + x2*z0 - x2*z1)/Delta;

  grad_phi_z[0] = -(x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2)/Delta;
  grad_phi_z[1] =  (x0*y2 - x0*y3 - x2*y0 + x2*y3 + x3*y0 - x3*y2)/Delta;
  grad_phi_z[2] = -(x0*y1 - x0*y3 - x1*y0 + x1*y3 + x3*y0 - x3*y1)/Delta;
  grad_phi_z[3] =  (x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)/Delta;

  return 0;
}

PetscErrorCode Calc3DP2_GradPhi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp,
  PetscScalar *grad_phi_x, PetscScalar *grad_phi_y, PetscScalar *grad_phi_z)
{
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];
  const PetscInt n3 = _inter[en][3];

  const PetscScalar x0 = _x[n0], y0 = _y[n0], z0 = _z[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1], z1 = _z[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2], z2 = _z[n2];
  const PetscScalar x3 = _x[n3], y3 = _y[n3], z3 = _z[n3];

  const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 + x1*y0*z2 
            -x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3  
            +x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0  
            +x3*y1*z2 + x3*y2*z0 - x3*y2*z1;
  const PetscScalar Delta2 = Delta*Delta;

  const PetscScalar C0yz =  ( y1*z2 - y1*z3 - y2*z1 + y2*z3 + y3*z1 - y3*z2);
  const PetscScalar C1yz =  (-y0*z2 + y0*z3 + y2*z0 - y2*z3 - y3*z0 + y3*z2);
  const PetscScalar C2yz = -(-y0*z1 + y0*z3 + y1*z0 - y1*z3 - y3*z0 + y3*z1);
  const PetscScalar C3yz =  (-y0*z1 + y0*z2 + y1*z0 - y1*z2 - y2*z0 + y2*z1);
  const PetscScalar C0xz = -(x1*z2 - x1*z3 - x2*z1 + x2*z3 + x3*z1 - x3*z2);
  const PetscScalar C1xz =  (x0*z2 - x0*z3 - x2*z0 + x2*z3 + x3*z0 - x3*z2);
  const PetscScalar C2xz = -(x0*z1 - x0*z3 - x1*z0 + x1*z3 + x3*z0 - x3*z1);
  const PetscScalar C3xz =  (x0*z1 - x0*z2 - x1*z0 + x1*z2 + x2*z0 - x2*z1);

  const PetscScalar C0xy =  (x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2);
  const PetscScalar C1xy = -(x0*y2 - x0*y3 - x2*y0 + x2*y3 + x3*y0 - x3*y2);
  const PetscScalar C2xy =  (x0*y1 - x0*y3 - x1*y0 + x1*y3 + x3*y0 - x3*y1);
  const PetscScalar C3xy = -(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1);

  const PetscScalar C0xyz = (-x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
  const PetscScalar C1xyz = ( x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0);
  const PetscScalar C2xyz = (-x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0);
  const PetscScalar C3xyz = ( x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0);

  grad_phi_x[0] = C0yz*(4*C0xy*zp + 4*C0xyz + 4*C0xz*yp + 4*C0yz*xp + Delta)/Delta2;
  grad_phi_x[1] = C1yz*(4*C1xy*zp + 4*C1xyz + 4*C1xz*yp + 4*C1yz*xp + Delta)/Delta2;
  grad_phi_x[2] = C2yz*(4*C2xy*zp + 4*C2xyz + 4*C2xz*yp + 4*C2yz*xp + Delta)/Delta2;
  grad_phi_x[3] = C3yz*(4*C3xy*zp + 4*C3xyz + 4*C3xz*yp + 4*C3yz*xp + Delta)/Delta2;
  grad_phi_x[4] = 4*(C0xy*C1yz*zp + C0xyz*C1yz + C0xz*C1yz*yp + C0yz*C1xy*zp + C0yz*C1xyz + C0yz*C1xz*yp + 2*C0yz*C1yz*xp)/Delta2;
  grad_phi_x[5] = 4*(C1xy*C2yz*zp + C1xyz*C2yz + C1xz*C2yz*yp + C1yz*C2xy*zp + C1yz*C2xyz + C1yz*C2xz*yp + 2*C1yz*C2yz*xp)/Delta2;
  grad_phi_x[6] = 4*(C0xy*C2yz*zp + C0xyz*C2yz + C0xz*C2yz*yp + C0yz*C2xy*zp + C0yz*C2xyz + C0yz*C2xz*yp + 2*C0yz*C2yz*xp)/Delta2;
  grad_phi_x[7] = 4*(C0xy*C3yz*zp + C0xyz*C3yz + C0xz*C3yz*yp + C0yz*C3xy*zp + C0yz*C3xyz + C0yz*C3xz*yp + 2*C0yz*C3yz*xp)/Delta2;
  grad_phi_x[8] = 4*(C1xy*C3yz*zp + C1xyz*C3yz + C1xz*C3yz*yp + C1yz*C3xy*zp + C1yz*C3xyz + C1yz*C3xz*yp + 2*C1yz*C3yz*xp)/Delta2;
  grad_phi_x[9] = 4*(C2xy*C3yz*zp + C2xyz*C3yz + C2xz*C3yz*yp + C2yz*C3xy*zp + C2yz*C3xyz + C2yz*C3xz*yp + 2*C2yz*C3yz*xp)/Delta2;

  grad_phi_y[0] = C0xz*(4*C0xy*zp + 4*C0xyz + 4*C0xz*yp + 4*C0yz*xp + Delta)/Delta2;
  grad_phi_y[1] = C1xz*(4*C1xy*zp + 4*C1xyz + 4*C1xz*yp + 4*C1yz*xp + Delta)/Delta2;
  grad_phi_y[2] = C2xz*(4*C2xy*zp + 4*C2xyz + 4*C2xz*yp + 4*C2yz*xp + Delta)/Delta2;
  grad_phi_y[3] = C3xz*(4*C3xy*zp + 4*C3xyz + 4*C3xz*yp + 4*C3yz*xp + Delta)/Delta2;
  grad_phi_y[4] = 4*(C0xy*C1xz*zp + C0xyz*C1xz + C0xz*C1xy*zp + C0xz*C1xyz + 2*C0xz*C1xz*yp + C0xz*C1yz*xp + C0yz*C1xz*xp)/Delta2;
  grad_phi_y[5] = 4*(C1xy*C2xz*zp + C1xyz*C2xz + C1xz*C2xy*zp + C1xz*C2xyz + 2*C1xz*C2xz*yp + C1xz*C2yz*xp + C1yz*C2xz*xp)/Delta2;
  grad_phi_y[6] = 4*(C0xy*C2xz*zp + C0xyz*C2xz + C0xz*C2xy*zp + C0xz*C2xyz + 2*C0xz*C2xz*yp + C0xz*C2yz*xp + C0yz*C2xz*xp)/Delta2;
  grad_phi_y[7] = 4*(C0xy*C3xz*zp + C0xyz*C3xz + C0xz*C3xy*zp + C0xz*C3xyz + 2*C0xz*C3xz*yp + C0xz*C3yz*xp + C0yz*C3xz*xp)/Delta2;
  grad_phi_y[8] = 4*(C1xy*C3xz*zp + C1xyz*C3xz + C1xz*C3xy*zp + C1xz*C3xyz + 2*C1xz*C3xz*yp + C1xz*C3yz*xp + C1yz*C3xz*xp)/Delta2;
  grad_phi_y[9] = 4*(C2xy*C3xz*zp + C2xyz*C3xz + C2xz*C3xy*zp + C2xz*C3xyz + 2*C2xz*C3xz*yp + C2xz*C3yz*xp + C2yz*C3xz*xp)/Delta2;

  grad_phi_z[0] = C0xy*(4*C0xy*zp + 4*C0xyz + 4*C0xz*yp + 4*C0yz*xp + Delta)/Delta2;
  grad_phi_z[1] = C1xy*(4*C1xy*zp + 4*C1xyz + 4*C1xz*yp + 4*C1yz*xp + Delta)/Delta2;
  grad_phi_z[2] = C2xy*(4*C2xy*zp + 4*C2xyz + 4*C2xz*yp + 4*C2yz*xp + Delta)/Delta2;
  grad_phi_z[3] = C3xy*(4*C3xy*zp + 4*C3xyz + 4*C3xz*yp + 4*C3yz*xp + Delta)/Delta2;
  grad_phi_z[4] = 4*(2*C0xy*C1xy*zp + C0xy*C1xyz + C0xy*C1xz*yp + C0xy*C1yz*xp + C0xyz*C1xy + C0xz*C1xy*yp + C0yz*C1xy*xp)/Delta2;
  grad_phi_z[5] = 4*(2*C1xy*C2xy*zp + C1xy*C2xyz + C1xy*C2xz*yp + C1xy*C2yz*xp + C1xyz*C2xy + C1xz*C2xy*yp + C1yz*C2xy*xp)/Delta2;
  grad_phi_z[6] = 4*(2*C0xy*C2xy*zp + C0xy*C2xyz + C0xy*C2xz*yp + C0xy*C2yz*xp + C0xyz*C2xy + C0xz*C2xy*yp + C0yz*C2xy*xp)/Delta2;
  grad_phi_z[7] = 4*(2*C0xy*C3xy*zp + C0xy*C3xyz + C0xy*C3xz*yp + C0xy*C3yz*xp + C0xyz*C3xy + C0xz*C3xy*yp + C0yz*C3xy*xp)/Delta2;
  grad_phi_z[8] = 4*(2*C1xy*C3xy*zp + C1xy*C3xyz + C1xy*C3xz*yp + C1xy*C3yz*xp + C1xyz*C3xy + C1xz*C3xy*yp + C1yz*C3xy*xp)/Delta2;
  grad_phi_z[9] = 4*(2*C2xy*C3xy*zp + C2xy*C3xyz + C2xy*C3xz*yp + C2xy*C3yz*xp + C2xyz*C3xy + C2xz*C3xy*yp + C2yz*C3xy*xp)/Delta2;

  return 0;
}

PetscErrorCode Calc3DP1b_GradPhi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp,
  PetscScalar *grad_phi_x, PetscScalar *grad_phi_y, PetscScalar *grad_phi_z)
{
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];
  const PetscInt n3 = _inter[en][3];

  const PetscScalar x0 = _x[n0], y0 = _y[n0], z0 = _z[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1], z1 = _z[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2], z2 = _z[n2];
  const PetscScalar x3 = _x[n3], y3 = _y[n3], z3 = _z[n3];

  const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 + x1*y0*z2 
            -x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3  
            +x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0  
            +x3*y1*z2 + x3*y2*z0 - x3*y2*z1;
  const PetscScalar Delta4 = Delta*Delta*Delta*Delta;

  const PetscScalar C0yz =  ( y1*z2 - y1*z3 - y2*z1 + y2*z3 + y3*z1 - y3*z2);
  const PetscScalar C1yz =  (-y0*z2 + y0*z3 + y2*z0 - y2*z3 - y3*z0 + y3*z2);
  const PetscScalar C2yz = -(-y0*z1 + y0*z3 + y1*z0 - y1*z3 - y3*z0 + y3*z1);
  const PetscScalar C3yz =  (-y0*z1 + y0*z2 + y1*z0 - y1*z2 - y2*z0 + y2*z1);
  const PetscScalar C0xz = -(x1*z2 - x1*z3 - x2*z1 + x2*z3 + x3*z1 - x3*z2);
  const PetscScalar C1xz =  (x0*z2 - x0*z3 - x2*z0 + x2*z3 + x3*z0 - x3*z2);
  const PetscScalar C2xz = -(x0*z1 - x0*z3 - x1*z0 + x1*z3 + x3*z0 - x3*z1);
  const PetscScalar C3xz =  (x0*z1 - x0*z2 - x1*z0 + x1*z2 + x2*z0 - x2*z1);

  const PetscScalar C0xy =  (x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2);
  const PetscScalar C1xy = -(x0*y2 - x0*y3 - x2*y0 + x2*y3 + x3*y0 - x3*y2);
  const PetscScalar C2xy =  (x0*y1 - x0*y3 - x1*y0 + x1*y3 + x3*y0 - x3*y1);
  const PetscScalar C3xy = -(x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1);

  const PetscScalar C0xyz = (-x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1);
  const PetscScalar C1xyz = ( x0*y2*z3 - x0*y3*z2 - x2*y0*z3 + x2*y3*z0 + x3*y0*z2 - x3*y2*z0);
  const PetscScalar C2xyz = (-x0*y1*z3 + x0*y3*z1 + x1*y0*z3 - x1*y3*z0 - x3*y0*z1 + x3*y1*z0);
  const PetscScalar C3xyz = ( x0*y1*z2 - x0*y2*z1 - x1*y0*z2 + x1*y2*z0 + x2*y0*z1 - x2*y1*z0);

  grad_phi_x[0] =  -C0yz/Delta + 64*C0yz*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C1yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 64*C2yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C3yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                                         
  grad_phi_x[1] =  64*C0yz*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - C1yz/Delta + 64*C1yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 64*C2yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C3yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                                           
  grad_phi_x[2] =  64*C0yz*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C1yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - C2yz/Delta - 64*C2yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C3yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                           
  grad_phi_x[3] =  64*C0yz*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C1yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 64*C2yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - C3yz/Delta + 64*C3yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                            
  grad_phi_x[4] =  -256*C0yz*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 256*C1yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 256*C2yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 256*C3yz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;

  grad_phi_y[0] =  -C0xz/Delta + 64*C0xz*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C1xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 64*C2xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C3xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                                      
  grad_phi_y[1] =  64*C0xz*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - C1xz/Delta + 64*C1xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 64*C2xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C3xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                                           
  grad_phi_y[2] =  64*C0xz*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C1xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - C2xz/Delta - 64*C2xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C3xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                          
  grad_phi_y[3] =  64*C0xz*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C1xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 64*C2xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - C3xz/Delta + 64*C3xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                            
  grad_phi_y[4] =  -256*C0xz*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 256*C1xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 256*C2xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 256*C3xz*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;

  grad_phi_z[0] =  -C0xy/Delta + 64*C0xy*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C1xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 64*C2xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C3xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                                              
  grad_phi_z[1] =  64*C0xy*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - C1xy/Delta + 64*C1xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 64*C2xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C3xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                                              
  grad_phi_z[2] =  64*C0xy*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C1xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - C2xy/Delta - 64*C2xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C3xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                           
  grad_phi_z[3] =  64*C0xy*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 64*C1xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 64*C2xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - C3xy/Delta + 64*C3xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ;                                                            
  grad_phi_z[4] =  -256*C0xy*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 256*C1xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 + 256*C2xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4 - 256*C3xy*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta4 ; 

  return 0;
}
#endif