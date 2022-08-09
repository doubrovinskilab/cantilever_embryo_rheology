#ifndef ElementalBases_C
#define ElementalBases_C

PetscErrorCode Calc2DP1_Phi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *phi)
{
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];

  const PetscScalar x0 = _x[n0], y0 = _y[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2];

  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;

  phi[0] =  ( xp*y1 - xp*y2 - x1*yp + x1*y2 + x2*yp - x2*y1)/Delta;
  phi[1] =  (-xp*y0 + xp*y2 + x0*yp - x0*y2 - x2*yp + x2*y0)/Delta;
  phi[2] = -(-xp*y0 + xp*y1 + x0*yp - x0*y1 - x1*yp + x1*y0)/Delta;
  
  return 0;
}

PetscErrorCode Calc2DP1b_Phi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *phi)
{
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];

  const PetscScalar x0 = _x[n0], y0 = _y[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2];

  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;
  const PetscScalar Delta3 = Delta*Delta*Delta;

  phi[0] =  (xp*y1 - xp*y2 - x1*yp + x1*y2 + x2*yp - x2*y1)/Delta + 9*(-xp*y0 + xp*y1 + x0*yp - x0*y1 - x1*yp + x1*y0)*(-xp*y0 + xp*y2 + x0*yp - x0*y2 - x2*yp + x2*y0)*(xp*y1 - xp*y2 - x1*yp + x1*y2 + x2*yp - x2*y1)/Delta3;
  phi[1] =  (-xp*y0 + xp*y2 + x0*yp - x0*y2 - x2*yp + x2*y0)/Delta + 9*(-xp*y0 + xp*y1 + x0*yp - x0*y1 - x1*yp + x1*y0)*(-xp*y0 + xp*y2 + x0*yp - x0*y2 - x2*yp + x2*y0)*(xp*y1 - xp*y2 - x1*yp + x1*y2 + x2*yp - x2*y1)/Delta3;
  phi[2] = -(-xp*y0 + xp*y1 + x0*yp - x0*y1 - x1*yp + x1*y0)/Delta + 9*(-xp*y0 + xp*y1 + x0*yp - x0*y1 - x1*yp + x1*y0)*(-xp*y0 + xp*y2 + x0*yp - x0*y2 - x2*yp + x2*y0)*(xp*y1 - xp*y2 - x1*yp + x1*y2 + x2*yp - x2*y1)/Delta3;
  phi[3] = -27*(-xp*y0 + xp*y1 + x0*yp - x0*y1 - x1*yp + x1*y0)*(-xp*y0 + xp*y2 + x0*yp - x0*y2 - x2*yp + x2*y0)*(xp*y1 - xp*y2 - x1*yp + x1*y2 + x2*yp - x2*y1)/Delta3;
  
  return 0;
}

PetscErrorCode Calc2DP2_Phi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *phi)
{   
  const PetscInt n0  = _inter[en][0];
  const PetscInt n1  = _inter[en][1];
  const PetscInt n2  = _inter[en][2];

  const PetscScalar x0 = _x[n0], y0 = _y[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2];

  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;
  const PetscScalar Delta2 = Delta*Delta;

  phi[0] = (    ( xp*y1 - xp*y2 - x1*yp + x1*y2 + x2*yp - x2*y1)*(-Delta + 2.*xp*y1 - 2.*xp*y2 - 2.*x1*yp + 2.*x1*y2 + 2.*x2*yp - 2.*x2*y1)/(Delta2));
  phi[1] = (    (-xp*y0 + xp*y2 + x0*yp - x0*y2 - x2*yp + x2*y0)*(-Delta - 2.*xp*y0 + 2.*xp*y2 + 2.*x0*yp - 2.*x0*y2 - 2.*x2*yp + 2.*x2*y0)/(Delta2));
  phi[2] = (    (-xp*y0 + xp*y1 + x0*yp - x0*y1 - x1*yp + x1*y0)*( Delta - 2.*xp*y0 + 2.*xp*y1 + 2.*x0*yp - 2.*x0*y1 - 2.*x1*yp + 2.*x1*y0)/(Delta2));
  phi[3] = ( 4.*(-xp*y0 + xp*y2 + x0*yp - x0*y2 - x2*yp + x2*y0)*( xp*y1 - xp*y2 - x1*yp + x1*y2 + x2*yp - x2*y1)/(Delta2));
  phi[4] = (-4.*(-xp*y0 + xp*y1 + x0*yp - x0*y1 - x1*yp + x1*y0)*(-xp*y0 + xp*y2 + x0*yp - x0*y2 - x2*yp + x2*y0)/(Delta2));
  phi[5] = (-4.*(-xp*y0 + xp*y1 + x0*yp - x0*y1 - x1*yp + x1*y0)*( xp*y1 - xp*y2 - x1*yp + x1*y2 + x2*yp - x2*y1)/(Delta2));

  return 0;
}

PetscErrorCode Calc3DP1_Phi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *phi)
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
  
  phi[0] = -(xp*y1*z2 - xp*y1*z3 - xp*y2*z1 + xp*y2*z3 + xp*y3*z1 - xp*y3*z2 - x1*xp*z2 
    + x1*xp*z3 + x1*y2*zp- x1*y2*z3 - x1*y3*zp+ x1*y3*z2 + x2*xp*z1 - x2*xp*z3 - x2*y1*zp+ x2*y1*z3 + x2*y3*zp
    - x2*y3*z1 - x3*xp*z1 + x3*xp*z2 + x3*y1*zp- x3*y1*z2 - x3*y2*zp+ x3*y2*z1)/Delta;
  phi[1] = -(-xp*y0*z2 + xp*y0*z3 + xp*y2*z0 - xp*y2*z3 - xp*y3*z0 + xp*y3*z2 + x0*xp*z2 - x0*xp*z3 
     - x0*y2*zp+ x0*y2*z3 + x0*y3*zp- x0*y3*z2 - x2*xp*z0 + x2*xp*z3 + x2*y0*zp- x2*y0*z3 - x2*y3*zp
     + x2*y3*z0 + x3*xp*z0 - x3*xp*z2 - x3*y0*zp+ x3*y0*z2 + x3*y2*zp- x3*y2*z0)/Delta;
  phi[2] =  (-xp*y0*z1 + xp*y0*z3 + xp*y1*z0 - xp*y1*z3 - xp*y3*z0 + xp*y3*z1 + x0*xp*z1 
    - x0*xp*z3 - x0*y1*zp+ x0*y1*z3 + x0*y3*zp- x0*y3*z1 - x1*xp*z0 + x1*xp*z3 + x1*y0*zp- x1*y0*z3 - x1*y3*zp
    + x1*y3*z0 + x3*xp*z0 - x3*xp*z1 - x3*y0*zp+ x3*y0*z1 + x3*y1*zp- x3*y1*z0)/Delta;
  phi[3] = -(-xp*y0*z1 + xp*y0*z2 + xp*y1*z0 - xp*y1*z2 - xp*y2*z0 + xp*y2*z1 + x0*xp*z1 
    - x0*xp*z2 - x0*y1*zp+ x0*y1*z2 + x0*y2*zp- x0*y2*z1 - x1*xp*z0 + x1*xp*z2 + x1*y0*zp- x1*y0*z2 - x1*y2*zp
    + x1*y2*z0 + x2*xp*z0 - x2*xp*z1 - x2*y0*zp+ x2*y0*z1 + x2*y1*zp- x2*y1*z0)/Delta;

  return 0;
}

PetscErrorCode Calc3DP2_Phi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *phi)
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

  phi[0] =    (C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(2.*C0xy*zp + 2.*C0xyz + 2.*C0xz*yp + 2.*C0yz*xp + Delta)/Delta2 ;
  phi[1] =    (C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(2.*C1xy*zp + 2.*C1xyz + 2.*C1xz*yp + 2.*C1yz*xp + Delta)/Delta2 ;
  phi[2] =    (C2xy*zp + C2xyz + C2xz*yp + C2yz*xp)*(2.*C2xy*zp + 2.*C2xyz + 2.*C2xz*yp + 2.*C2yz*xp + Delta)/Delta2 ;
  phi[3] =    (C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)*(2.*C3xy*zp + 2.*C3xyz + 2.*C3xz*yp + 2.*C3yz*xp + Delta)/Delta2 ;
  phi[4] = 4.*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(   C1xy*zp +    C1xyz +    C1xz*yp +    C1yz*xp)/Delta2 ;
  phi[5] = 4.*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(   C2xy*zp +    C2xyz +    C2xz*yp +    C2yz*xp)/Delta2 ;
  phi[6] = 4.*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(   C2xy*zp +    C2xyz +    C2xz*yp +    C2yz*xp)/Delta2 ;
  phi[7] = 4.*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(   C3xy*zp +    C3xyz +    C3xz*yp +    C3yz*xp)/Delta2 ;
  phi[8] = 4.*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(   C3xy*zp +    C3xyz +    C3xz*yp +    C3yz*xp)/Delta2 ;
  phi[9] = 4.*(C2xy*zp + C2xyz + C2xz*yp + C2yz*xp)*(   C3xy*zp +    C3xyz +    C3xz*yp +    C3yz*xp)/Delta2 ;

  return 0;
}

PetscErrorCode Calc3DP1b_Phi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z,
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *phi)
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

  phi[0] = -(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)/Delta + 64*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4;
  phi[1] = -(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)/Delta + 64*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4;
  phi[2] =  (-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)/Delta + 64*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4;
  phi[3] = -(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta + 64*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4;
  phi[4] = -256*(C0xy*zp + C0xyz + C0xz*yp + C0yz*xp)*(C1xy*zp + C1xyz + C1xz*yp + C1yz*xp)*(-C2xy*zp - C2xyz - C2xz*yp - C2yz*xp)*(C3xy*zp + C3xyz + C3xz*yp + C3yz*xp)/Delta4;

  return 0;
}
#endif