#ifndef FEMFunctions_H
#define FEMFunctions_H

PetscErrorCode CheckBaseCompatibilityFor2DP1(const PetscMPIInt rank, const PetscInt en, 
  PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter)
{
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];

  const PetscScalar x0 = _x[n0], y0 = _y[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2];

  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;

  PetscInt check = 0, n;
  PetscScalar base_0, base_1, base_2, x, y;
  for (PetscInt i=0; i < 3; i++)
  {
    n = _inter[en][i];
    x = _x[n];
    y = _y[n];

    base_0 = (PetscScalar) PetscAbsReal( ( x*y1 - x*y2 - x1*y + x1*y2 + x2*y - x2*y1)/Delta);
    base_1 = (PetscScalar) PetscAbsReal( (-x*y0 + x*y2 + x0*y - x0*y2 - x2*y + x2*y0)/Delta);
    base_2 = (PetscScalar) PetscAbsReal(-(-x*y0 + x*y1 + x0*y - x0*y1 - x1*y + x1*y0)/Delta);

    if ( ( i == 0 && base_0 > (1 - EPS) && base_1 < EPS && base_2 < EPS) || 
         ( i == 1 && base_0 < EPS && base_1 > (1 - EPS) && base_2 < EPS) ||
         ( i == 2 && base_0 < EPS && base_1 < EPS && base_2 > (1 - EPS))   )
      check++;
    else
      {
        SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The P1 triangle nodes are not in agreement with the derived bases.\n");
      }
  }
  if (check != 3)
  {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Not all P1 triangle nodes passed the bases test.\n");
  }

  return 0;
}


PetscErrorCode CheckBaseCompatibilityFor2DP1b(const PetscMPIInt rank, const PetscInt en, 
  PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter)
{
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];

  const PetscScalar x0 = _x[n0], y0 = _y[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2];

  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;
  const PetscScalar Delta3 = Delta*Delta*Delta;

  PetscInt check = 0, n;
  PetscScalar base_0, base_1, base_2, base_3, x, y;
  for (PetscInt i=0; i < 4; i++)
  {
    n = _inter[en][i];
    x = _x[n];
    y = _y[n];

    base_0 = (PetscScalar) PetscAbsReal( (x*y1 - x*y2 - x1*y + x1*y2 + x2*y - x2*y1)/Delta + 9*(-x*y0 + x*y1 + x0*y - x0*y1 - x1*y + x1*y0)*(-x*y0 + x*y2 + x0*y - x0*y2 - x2*y + x2*y0)*(x*y1 - x*y2 - x1*y + x1*y2 + x2*y - x2*y1)/Delta3);
    base_1 = (PetscScalar) PetscAbsReal( (-x*y0 + x*y2 + x0*y - x0*y2 - x2*y + x2*y0)/Delta + 9*(-x*y0 + x*y1 + x0*y - x0*y1 - x1*y + x1*y0)*(-x*y0 + x*y2 + x0*y - x0*y2 - x2*y + x2*y0)*(x*y1 - x*y2 - x1*y + x1*y2 + x2*y - x2*y1)/Delta3);
    base_2 = (PetscScalar) PetscAbsReal(-(-x*y0 + x*y1 + x0*y - x0*y1 - x1*y + x1*y0)/Delta + 9*(-x*y0 + x*y1 + x0*y - x0*y1 - x1*y + x1*y0)*(-x*y0 + x*y2 + x0*y - x0*y2 - x2*y + x2*y0)*(x*y1 - x*y2 - x1*y + x1*y2 + x2*y - x2*y1)/Delta3);
    base_3 = (PetscScalar) PetscAbsReal( -27*(-x*y0 + x*y1 + x0*y - x0*y1 - x1*y + x1*y0)*(-x*y0 + x*y2 + x0*y - x0*y2 - x2*y + x2*y0)*(x*y1 - x*y2 - x1*y + x1*y2 + x2*y - x2*y1)/Delta3);

    if (( i == 0 && base_0 > (1 - EPS) && base_1 < EPS && base_2 < EPS && base_3 < EPS ) || 
        ( i == 1 && base_0 < EPS && base_1 > (1 - EPS) && base_2 < EPS && base_3 < EPS ) ||
        ( i == 2 && base_0 < EPS && base_1 < EPS && base_2 > (1 - EPS) && base_3 < EPS ) ||
        ( i == 3 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 > (1 - EPS) ))
      check++;
    else
      {
        SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The P1b triangle nodes are not in agreement with the derived bases.\n");
      }
  }
  if (check != 4)
  {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Not all P1b triangle nodes passed the bases test.\n");
  }

  return 0;
}

PetscErrorCode CheckBaseCompatibilityFor2DP2(const PetscMPIInt rank, const PetscInt en, 
  PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter)
{
  const PetscInt n0  = _inter[en][0];
  const PetscInt n1  = _inter[en][1];
  const PetscInt n2  = _inter[en][2];

  const PetscScalar x0 = _x[n0], y0 = _y[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2];

  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;
  const PetscScalar Delta2 = Delta*Delta;

  PetscInt check = 0, n;
  PetscScalar base_0, base_1, base_2, base_3, base_4, base_5, x, y;
  for (PetscInt i=0; i < 6; i++)
  {
    n = _inter[en][i];
    x = _x[n];
    y = _y[n];

    base_0 = (PetscScalar) PetscAbsReal(   ( x*y1 - x*y2 - x1*y + x1*y2 + x2*y - x2*y1)*(-Delta + 2.*x*y1 - 2.*x*y2 - 2.*x1*y + 2.*x1*y2 + 2.*x2*y - 2.*x2*y1)/(Delta2));
    base_1 = (PetscScalar) PetscAbsReal(   (-x*y0 + x*y2 + x0*y - x0*y2 - x2*y + x2*y0)*(-Delta - 2.*x*y0 + 2.*x*y2 + 2.*x0*y - 2.*x0*y2 - 2.*x2*y + 2.*x2*y0)/(Delta2));
    base_2 = (PetscScalar) PetscAbsReal(   (-x*y0 + x*y1 + x0*y - x0*y1 - x1*y + x1*y0)*( Delta - 2.*x*y0 + 2.*x*y1 + 2.*x0*y - 2.*x0*y1 - 2.*x1*y + 2.*x1*y0)/(Delta2));
    base_3 = (PetscScalar) PetscAbsReal( 4.*(-x*y0 + x*y2 + x0*y - x0*y2 - x2*y + x2*y0)*(x*y1 - x*y2 - x1*y + x1*y2 + x2*y - x2*y1)/(Delta2));
    base_4 = (PetscScalar) PetscAbsReal(-4*(-x*y0 + x*y1 + x0*y - x0*y1 - x1*y + x1*y0)*(-x*y0 + x*y2 + x0*y - x0*y2 - x2*y + x2*y0)/(Delta2));
    base_5 = (PetscScalar) PetscAbsReal(-4*(-x*y0 + x*y1 + x0*y - x0*y1 - x1*y + x1*y0)*(x*y1 - x*y2 - x1*y + x1*y2 + x2*y - x2*y1)/(Delta2));

    if (( i == 0 && base_0 > (1 - EPS) && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 < EPS && base_5 < EPS) || 
        ( i == 1 && base_0 < EPS && base_1 > (1 - EPS) && base_2 < EPS && base_3 < EPS && base_4 < EPS && base_5 < EPS) ||
        ( i == 2 && base_0 < EPS && base_1 < EPS && base_2 > (1 - EPS) && base_3 < EPS && base_4 < EPS && base_5 < EPS) ||
        ( i == 3 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 > (1 - EPS) && base_4 < EPS && base_5 < EPS) ||
        ( i == 4 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 > (1 - EPS) && base_5 < EPS) ||
        ( i == 5 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 < EPS && base_5 > (1 - EPS))  )
      check++;
    else
      {
        SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The P2 triangle nodes are not in agreement with the derived bases.\n");
      }
  }
  if (check != 6)
  {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Not all P2 triangle nodes passed the bases test.\n");
  }

  return 0;
}

PetscErrorCode CheckBaseCompatibilityFor3DP1(const PetscMPIInt rank, const PetscInt en, 
  PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter)
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

  PetscInt check = 0, n;
  PetscScalar base_0, base_1, base_2, base_3, x, y, z;
  for (PetscInt i=0; i < 4; i++)
  {
    n = _inter[en][i];
    x = _x[n], y = _y[n], z = _z[n];

    base_0 = (PetscScalar) PetscAbsReal(-(x*y1*z2 - x*y1*z3 - x*y2*z1 + x*y2*z3 + x*y3*z1 - x*y3*z2 - x1*y*z2 + x1*y*z3 
     + x1*y2*z - x1*y2*z3 - x1*y3*z + x1*y3*z2 + x2*y*z1 - x2*y*z3 - x2*y1*z + x2*y1*z3 + x2*y3*z 
     - x2*y3*z1 - x3*y*z1 + x3*y*z2 + x3*y1*z - x3*y1*z2 - x3*y2*z + x3*y2*z1)/Delta);

    base_1 = (PetscScalar) PetscAbsReal(-(-x*y0*z2 + x*y0*z3 + x*y2*z0 - x*y2*z3 - x*y3*z0 + x*y3*z2 + x0*y*z2 - x0*y*z3 
     - x0*y2*z + x0*y2*z3 + x0*y3*z - x0*y3*z2 - x2*y*z0 + x2*y*z3 + x2*y0*z - x2*y0*z3 - x2*y3*z 
     + x2*y3*z0 + x3*y*z0 - x3*y*z2 - x3*y0*z + x3*y0*z2 + x3*y2*z - x3*y2*z0)/Delta);

    base_2 = (PetscScalar) PetscAbsReal( (-x*y0*z1 + x*y0*z3 + x*y1*z0 - x*y1*z3 - x*y3*z0 + x*y3*z1 + x0*y*z1 - x0*y*z3 
     - x0*y1*z + x0*y1*z3 + x0*y3*z - x0*y3*z1 - x1*y*z0 + x1*y*z3 + x1*y0*z - x1*y0*z3 - x1*y3*z 
     + x1*y3*z0 + x3*y*z0 - x3*y*z1 - x3*y0*z + x3*y0*z1 + x3*y1*z - x3*y1*z0)/Delta);

    base_3 = (PetscScalar) PetscAbsReal(-(-x*y0*z1 + x*y0*z2 + x*y1*z0 - x*y1*z2 - x*y2*z0 + x*y2*z1 + x0*y*z1 - x0*y*z2 
     - x0*y1*z + x0*y1*z2 + x0*y2*z - x0*y2*z1 - x1*y*z0 + x1*y*z2 + x1*y0*z - x1*y0*z2 - x1*y2*z 
     + x1*y2*z0 + x2*y*z0 - x2*y*z1 - x2*y0*z + x2*y0*z1 + x2*y1*z - x2*y1*z0)/Delta);

    if ( (i == 0 && base_0 > (1 - EPS) && base_1 < EPS && base_2 < EPS && base_3 < EPS) || 
         (i == 1 && base_0 < EPS && base_1 > (1 - EPS) && base_2 < EPS && base_3 < EPS) ||
         (i == 2 && base_0 < EPS && base_1 < EPS && base_2 > (1 - EPS) && base_3 < EPS) ||
         (i == 3 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 > (1 - EPS)) )
      check++;
    else
      {
        SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The P1 tetraherdal nodes are not in agreement with the derived bases.\n");
      }
  }

  if (check != 4)
  {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Not all P1 tetraherdal nodes passed the bases test.\n");
  }
    
  return 0;
}


PetscErrorCode CheckBaseCompatibilityFor3DP1b(const PetscMPIInt rank, const PetscInt en, 
  PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter)
{
  const PetscInt n0 = _inter[en][0];
  const PetscInt n1 = _inter[en][1];
  const PetscInt n2 = _inter[en][2];
  const PetscInt n3 = _inter[en][3];

  const PetscScalar x0 = _x[n0], y0 = _y[n0], z0 = _z[n0];
  const PetscScalar x1 = _x[n1], y1 = _y[n1], z1 = _z[n1];
  const PetscScalar x2 = _x[n2], y2 = _y[n2], z2 = _z[n2];
  const PetscScalar x3 = _x[n3], y3 = _y[n3], z3 = _z[n3];

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

  const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2  
    + x1*y0*z2 - x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3 
    + x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0 + x3*y1*z2 
    + x3*y2*z0 - x3*y2*z1;
  const PetscScalar Delta4 = Delta*Delta*Delta*Delta;

  PetscInt check = 0, n;
  PetscScalar base_0, base_1, base_2, base_3, base_4, x, y, z;
  for (PetscInt i=0; i < 5; i++)
  {
    n = _inter[en][i];
    x = _x[n], y = _y[n], z = _z[n];

    base_0 =  (PetscScalar) PetscAbsReal(-(C0xy*z + C0xyz + C0xz*y + C0yz*x)/Delta + 
      64.*(C0xy*z + C0xyz + C0xz*y + C0yz*x)*(C1xy*z + C1xyz + C1xz*y + C1yz*x)*(-C2xy*z - C2xyz - C2xz*y - C2yz*x)*(C3xy*z + C3xyz + C3xz*y + C3yz*x)/Delta4);
    base_1 =  (PetscScalar) PetscAbsReal(-(C1xy*z + C1xyz + C1xz*y + C1yz*x)/Delta + 
      64.*(C0xy*z + C0xyz + C0xz*y + C0yz*x)*(C1xy*z + C1xyz + C1xz*y + C1yz*x)*(-C2xy*z - C2xyz - C2xz*y - C2yz*x)*(C3xy*z + C3xyz + C3xz*y + C3yz*x)/Delta4);
    base_2 =  (PetscScalar) PetscAbsReal((-C2xy*z - C2xyz - C2xz*y - C2yz*x)/Delta + 
      64.*(C0xy*z + C0xyz + C0xz*y + C0yz*x)*(C1xy*z + C1xyz + C1xz*y + C1yz*x)*(-C2xy*z - C2xyz - C2xz*y - C2yz*x)*(C3xy*z + C3xyz + C3xz*y + C3yz*x)/Delta4);
    base_3 =  (PetscScalar) PetscAbsReal(-(C3xy*z + C3xyz + C3xz*y + C3yz*x)/Delta + 
      64.*(C0xy*z + C0xyz + C0xz*y + C0yz*x)*(C1xy*z + C1xyz + C1xz*y + C1yz*x)*(-C2xy*z - C2xyz - C2xz*y - C2yz*x)*(C3xy*z + C3xyz + C3xz*y + C3yz*x)/Delta4);
    base_4 =  (PetscScalar) PetscAbsReal(-256.*(C0xy*z + C0xyz + C0xz*y + C0yz*x)*(C1xy*z + C1xyz + C1xz*y + C1yz*x)*(-C2xy*z - C2xyz - C2xz*y - C2yz*x)*(C3xy*z + C3xyz + C3xz*y + C3yz*x)/Delta4);          

    if ( (i == 0 && base_0 > (1 - EPS) && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 < EPS) || 
         (i == 1 && base_0 < EPS && base_1 > (1 - EPS) && base_2 < EPS && base_3 < EPS && base_4 < EPS) ||
         (i == 2 && base_0 < EPS && base_1 < EPS && base_2 > (1 - EPS) && base_3 < EPS && base_4 < EPS) ||
         (i == 3 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 > (1 - EPS) && base_4 < EPS) ||
         (i == 4 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 > (1 - EPS) )   )
      check++;
    else
      {
        SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The P1b tetraherdal nodes are not in agreement with the derived bases.\n");
      }
  }

  if (check != 5)
  {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Not all P1 tetraherdal nodes passed the bases test.\n");
  }
    
  return 0;
}

PetscErrorCode CheckBaseCompatibilityFor3DP2(const PetscMPIInt rank, const PetscInt en, 
  PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter)
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

  PetscInt check = 0, n;
  PetscScalar base_0, base_1, base_2, base_3, base_4, base_5, base_6, base_7, base_8, base_9, x, y, z;
  for (PetscInt i=0; i < 10; i++)
  {
    n = _inter[en][i];
    x = _x[n], y = _y[n], z = _z[n];

    base_0 = (PetscScalar) PetscAbsReal(    (C0xy*z + C0xyz + C0xz*y + C0yz*x)*(2*C0xy*z + 2.*C0xyz + 2.*C0xz*y + 2.*C0yz*x + Delta)/Delta2 );
    base_1 = (PetscScalar) PetscAbsReal(    (C1xy*z + C1xyz + C1xz*y + C1yz*x)*(2*C1xy*z + 2.*C1xyz + 2.*C1xz*y + 2.*C1yz*x + Delta)/Delta2 );
    base_2 = (PetscScalar) PetscAbsReal(    (C2xy*z + C2xyz + C2xz*y + C2yz*x)*(2*C2xy*z + 2.*C2xyz + 2.*C2xz*y + 2.*C2yz*x + Delta)/Delta2 );
    base_3 = (PetscScalar) PetscAbsReal(    (C3xy*z + C3xyz + C3xz*y + C3yz*x)*(2*C3xy*z + 2.*C3xyz + 2.*C3xz*y + 2.*C3yz*x + Delta)/Delta2 );
    base_4 = (PetscScalar) PetscAbsReal( 4.*(C0xy*z + C0xyz + C0xz*y + C0yz*x)*(C1xy*z + C1xyz + C1xz*y + C1yz*x)/Delta2 );
    base_5 = (PetscScalar) PetscAbsReal( 4.*(C1xy*z + C1xyz + C1xz*y + C1yz*x)*(C2xy*z + C2xyz + C2xz*y + C2yz*x)/Delta2 );
    base_6 = (PetscScalar) PetscAbsReal( 4.*(C0xy*z + C0xyz + C0xz*y + C0yz*x)*(C2xy*z + C2xyz + C2xz*y + C2yz*x)/Delta2 );
    base_7 = (PetscScalar) PetscAbsReal( 4.*(C0xy*z + C0xyz + C0xz*y + C0yz*x)*(C3xy*z + C3xyz + C3xz*y + C3yz*x)/Delta2 );
    base_8 = (PetscScalar) PetscAbsReal( 4.*(C1xy*z + C1xyz + C1xz*y + C1yz*x)*(C3xy*z + C3xyz + C3xz*y + C3yz*x)/Delta2 );
    base_9 = (PetscScalar) PetscAbsReal( 4.*(C2xy*z + C2xyz + C2xz*y + C2yz*x)*(C3xy*z + C3xyz + C3xz*y + C3yz*x)/Delta2 );

    if (( i == 0 && base_0 > (1 - EPS) && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 < EPS && base_5 < EPS && base_6 < EPS && base_7 < EPS && base_8 < EPS && base_9 < EPS) || 
        ( i == 1 && base_0 < EPS && base_1 > (1 - EPS) && base_2 < EPS && base_3 < EPS && base_4 < EPS && base_5 < EPS && base_6 < EPS && base_7 < EPS && base_8 < EPS && base_9 < EPS) ||
        ( i == 2 && base_0 < EPS && base_1 < EPS && base_2 > (1 - EPS) && base_3 < EPS && base_4 < EPS && base_5 < EPS && base_6 < EPS && base_7 < EPS && base_8 < EPS && base_9 < EPS) ||
        ( i == 3 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 > (1 - EPS) && base_4 < EPS && base_5 < EPS && base_6 < EPS && base_7 < EPS && base_8 < EPS && base_9 < EPS) ||
        ( i == 4 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 > (1 - EPS) && base_5 < EPS && base_6 < EPS && base_7 < EPS && base_8 < EPS && base_9 < EPS) ||
        ( i == 5 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 < EPS && base_5 > (1 - EPS) && base_6 < EPS && base_7 < EPS && base_8 < EPS && base_9 < EPS) ||
        ( i == 6 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 < EPS && base_5 < EPS && base_6 > (1 - EPS) && base_7 < EPS && base_8 < EPS && base_9 < EPS) ||
        ( i == 7 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 < EPS && base_5 < EPS && base_6 < EPS && base_7 > (1 - EPS) && base_8 < EPS && base_9 < EPS) ||
        ( i == 8 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 < EPS && base_5 < EPS && base_6 < EPS && base_7 < EPS && base_8 > (1 - EPS) && base_9 < EPS) ||
        ( i == 9 && base_0 < EPS && base_1 < EPS && base_2 < EPS && base_3 < EPS && base_4 < EPS && base_5 < EPS && base_6 < EPS && base_7 < EPS && base_8 < EPS && base_9 > (1 - EPS) ))
      check++;
    else
      {
        SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The P2 tetraherdal nodes are not in agreement with the derived bases.\n");
      }
  }

  if (check != 10)
  {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Not all P2 tetraherdal nodes passed the bases test.\n");
  }
  return 0;
}

PetscErrorCode CheckBaseCompatibility(System_Struct *system, FLUID_ELEMENT_TYPE Elem_type,
  const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter)
{
  if (system->Dim == 2)
  {
    if (Elem_type == P1)
      CheckBaseCompatibilityFor2DP1(system->rank, en, _x, _y, _z, _inter);
    else if (Elem_type == P1b)
      CheckBaseCompatibilityFor2DP1b(system->rank, en, _x, _y, _z, _inter);
    else if (Elem_type == P2)
      CheckBaseCompatibilityFor2DP2(system->rank,  en, _x, _y, _z, _inter);
  }
  else if (system->Dim == 3)
  {
    if (Elem_type == P1)
      CheckBaseCompatibilityFor3DP1(system->rank, en, _x, _y, _z, _inter);
    else if (Elem_type == P1b)
      CheckBaseCompatibilityFor3DP1b(system->rank, en, _x, _y, _z, _inter);
    else if (Elem_type == P2)
      CheckBaseCompatibilityFor3DP2(system->rank, en, _x, _y, _z, _inter);
  }
  else
  {
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! Incorrect dimensions");          
  } 

  return 0;
}


PetscErrorCode CheckBoundaryNodesArrangementFor2DP2(const PetscMPIInt rank, const PetscInt bn, 
  PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_bound)
{
  const PetscInt n0   = _bound[bn][0];
  const PetscInt n1   = _bound[bn][1];
  const PetscInt n01  = _bound[bn][2];

  const PetscScalar x0  = _x[n0], y0 = _y[n0];
  const PetscScalar x1  = _x[n1], y1 = _y[n1];
  const PetscScalar x01 = _x[n01], y01 = _y[n01];

  if ((PetscAbsReal(x0+x1 - 2*x01) > EPS) || (PetscAbsReal(y0+y1 - 2*y01) > EPS)) {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The quadratic triangle boundary nodes are not in agreement with the derived bases.\n");
  }
  return 0;
}


PetscErrorCode CheckBoundaryNodesArrangementFor3DP2(const PetscMPIInt rank, const PetscInt bn, 
  PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_bound)
{
  const PetscInt n0  = _bound[bn][0];
  const PetscInt n1  = _bound[bn][1];
  const PetscInt n2  = _bound[bn][2];
  const PetscInt n01 = _bound[bn][3];
  const PetscInt n12 = _bound[bn][4];
  const PetscInt n02 = _bound[bn][5];

  const PetscScalar x0  = _x[n0], y0 = _y[n0], z0 = _z[n0];
  const PetscScalar x1  = _x[n1], y1 = _y[n1], z1 = _z[n1];
  const PetscScalar x2  = _x[n2], y2 = _y[n2], z2 = _z[n2];
  const PetscScalar x01 = _x[n01], y01 = _y[n01], z01 = _z[n01];
  const PetscScalar x12 = _x[n12], y12 = _y[n12], z12 = _z[n12];
  const PetscScalar x02 = _x[n02], y02 = _y[n02], z02 = _z[n02];

  if ((PetscAbsReal(x0+x1 - 2*x01) > EPS) || (PetscAbsReal(y0+y1 - 2*y01) > EPS) || (PetscAbsReal(z0+z1 - 2*z01) > EPS)) {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The quadratic tetrahedron boundary nodes are not in agreement with the derived bases.\n");
  }
  if ((PetscAbsReal(x1+x2 - 2*x12) > EPS) || (PetscAbsReal(y1+y2 - 2*y12) > EPS) || (PetscAbsReal(z1+z2 - 2*z12) > EPS)) {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The quadratic tetrahedron boundary nodes are not in agreement with the derived bases.\n");
  }
  if ((PetscAbsReal(x0+x2 - 2*x02) > EPS) || (PetscAbsReal(y0+y2 - 2*y02) > EPS) || (PetscAbsReal(z0+z2 - 2*z02) > EPS)) {
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"The quadratic tetrahedron boundary nodes are not in agreement with the derived bases.\n");
  }
  return 0;
}

PetscErrorCode CheckBoundaryCompatibility(System_Struct *system, FLUID_ELEMENT_TYPE Elem_type,
  const PetscInt N_bound_loc, const PetscInt bn, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_bound)
{
  if (system->Dim == 2)
  {
    if (Elem_type == P2)
      if (N_bound_loc > 0) 
        CheckBoundaryNodesArrangementFor2DP2(system->rank, bn, _x, _y, _z, _bound);
  }
  else if (system->Dim == 3)
  {
    if (Elem_type == P2)
      if (N_bound_loc > 0) 
        CheckBaseCompatibilityFor3DP2(system->rank, bn, _x, _y, _z, _bound);
  }
  else
  {
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! Incorrect dimensions");          
  } 

  return 0;
}

#include "ElementalFunctions.c"

#endif