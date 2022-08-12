#ifndef ElementalMatrixMass_C
#define ElementalMatrixMass_C

PetscErrorCode Calc2DP1_Mek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Mek)
{
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];
  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;
  Mek[0][0] =  Delta/12.;
  Mek[0][1] =  Delta/24.;
  Mek[0][2] =  Delta/24.;
  Mek[1][0] =  Delta/24.;
  Mek[1][1] =  Delta/12.;
  Mek[1][2] =  Delta/24.;
  Mek[2][0] =  Delta/24.;
  Mek[2][1] =  Delta/24.;
  Mek[2][2] =  Delta/12.;  

  return 0;
}


PetscErrorCode Calc2DP1b_Mek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Mek)
{
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];
  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;
  Mek[0][0] =  83.*Delta/1680.;
  Mek[0][1] =  13.*Delta/1680.;
  Mek[0][2] =  13.*Delta/1680.;
  Mek[0][3] =   3.*Delta/112.;
  Mek[1][0] =  13.*Delta/1680.;
  Mek[1][1] =  83.*Delta/1680.;
  Mek[1][2] =  13.*Delta/1680.;
  Mek[1][3] =   3.*Delta/112.;
  Mek[2][0] =  13.*Delta/1680.;
  Mek[2][1] =  13.*Delta/1680.;
  Mek[2][2] =  83.*Delta/1680.;
  Mek[2][3] =   3.*Delta/112.;
  Mek[3][0] =   3.*Delta/112.;
  Mek[3][1] =   3.*Delta/112.;
  Mek[3][2] =   3.*Delta/112.;
  Mek[3][3] =  81.*Delta/560.;
  return 0;
}

PetscErrorCode Calc2DP2_Mek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Mek)
{
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];

  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;

  Mek[0][0] =    Delta/60.;
  Mek[0][1] =  -Delta/360.;
  Mek[0][2] =  -Delta/360.;
  Mek[0][3] =           0.;
  Mek[0][4] =   -Delta/90.;
  Mek[0][5] =           0.;
  
  Mek[1][0] =  -Delta/360.;
  Mek[1][1] =    Delta/60.;
  Mek[1][2] =  -Delta/360.;
  Mek[1][3] =           0.;
  Mek[1][4] =           0.;
  Mek[1][5] =   -Delta/90.;

  Mek[2][0] =  -Delta/360.;
  Mek[2][1] =  -Delta/360.;
  Mek[2][2] =    Delta/60.;
  Mek[2][3] =   -Delta/90.;
  Mek[2][4] =           0.;
  Mek[2][5] =           0.;

  Mek[3][0] =           0.;
  Mek[3][1] =           0.;
  Mek[3][2] =   -Delta/90.;
  Mek[3][3] =  4*Delta/45.;
  Mek[3][4] =  2*Delta/45.;
  Mek[3][5] =  2*Delta/45.;

  Mek[4][0] =   -Delta/90.;
  Mek[4][1] =           0.;
  Mek[4][2] =           0.;
  Mek[4][3] =  2*Delta/45.;
  Mek[4][4] =  4*Delta/45.;
  Mek[4][5] =  2*Delta/45.;

  Mek[5][0] =   0.;
  Mek[5][1] =   -Delta/90.;
  Mek[5][2] =   0.;
  Mek[5][3] =  2*Delta/45.;
  Mek[5][4] =  2*Delta/45.;
  Mek[5][5] =  4*Delta/45.;

  return 0;
}


PetscErrorCode Calc3DP1_Mek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Mek)
{
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscInt n3 = e_inter[en][3];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];
  
  const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 + x1*y0*z2 
            -x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3  
            +x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0  
            +x3*y1*z2 + x3*y2*z0 - x3*y2*z1;

  Mek[0][0] =  Delta/60.;
  Mek[0][1] =  Delta/120.;
  Mek[0][2] =  Delta/120.;
  Mek[0][3] =  Delta/120.;
  Mek[1][0] =  Delta/120.;
  Mek[1][1] =  Delta/60.;
  Mek[1][2] =  Delta/120.;
  Mek[1][3] =  Delta/120.;
  Mek[2][0] =  Delta/120.;
  Mek[2][1] =  Delta/120.;
  Mek[2][2] =  Delta/60.;
  Mek[2][3] =  Delta/120.;
  Mek[3][0] =  Delta/120.;
  Mek[3][1] =  Delta/120.;
  Mek[3][2] =  Delta/120.;
  Mek[3][3] =  Delta/60.;

  return 0;
}


PetscErrorCode Calc3DP2_Mek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Mek)
{
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscInt n3 = e_inter[en][3];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];
  
  const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 + x1*y0*z2 
            -x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3  
            +x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0  
            +x3*y1*z2 + x3*y2*z0 - x3*y2*z1;

  Mek[0][0] =  Delta/420.;
  Mek[0][1] =  Delta/2520.;
  Mek[0][2] =  Delta/2520.;
  Mek[0][3] =  Delta/2520.;
  Mek[0][4] =  -Delta/630.;
  Mek[0][5] =  -Delta/420.;
  Mek[0][6] =  -Delta/630.;
  Mek[0][7] =  -Delta/630.;
  Mek[0][8] =  -Delta/420.;
  Mek[0][9] =  -Delta/420.;
  Mek[1][0] =  Delta/2520.;
  Mek[1][1] =  Delta/420.; 
  Mek[1][2] =  Delta/2520.;
  Mek[1][3] =  Delta/2520.;
  Mek[1][4] =  -Delta/630.;
  Mek[1][5] =  -Delta/630.;
  Mek[1][6] =  -Delta/420.;
  Mek[1][7] =  -Delta/420.;
  Mek[1][8] =  -Delta/630.;
  Mek[1][9] =  -Delta/420.;
  Mek[2][0] =  Delta/2520.;
  Mek[2][1] =  Delta/2520.;
  Mek[2][2] =  Delta/420.; 
  Mek[2][3] =  Delta/2520.;
  Mek[2][4] =  -Delta/420.;
  Mek[2][5] =  -Delta/630.;
  Mek[2][6] =  -Delta/630.;
  Mek[2][7] =  -Delta/420.;
  Mek[2][8] =  -Delta/420.;
  Mek[2][9] =  -Delta/630.;
  Mek[3][0] =  Delta/2520.;
  Mek[3][1] =  Delta/2520.;
  Mek[3][2] =  Delta/2520.;
  Mek[3][3] =  Delta/420.; 
  Mek[3][4] =  -Delta/420.;
  Mek[3][5] =  -Delta/420.;
  Mek[3][6] =  -Delta/420.;
  Mek[3][7] =  -Delta/630.;
  Mek[3][8] =  -Delta/630.;
  Mek[3][9] =  -Delta/630.;
  Mek[4][0] =  -Delta/630.;
  Mek[4][1] =  -Delta/630.;
  Mek[4][2] =  -Delta/420.;
  Mek[4][3] =  -Delta/420.; 
  Mek[4][4] =  4*Delta/315.;
  Mek[4][5] =  2*Delta/315.;
  Mek[4][6] =  2*Delta/315.;
  Mek[4][7] =  2*Delta/315.;
  Mek[4][8] =  2*Delta/315.;
  Mek[4][9] =  Delta/315.;
  Mek[5][0] =  -Delta/420.;
  Mek[5][1] =  -Delta/630.;
  Mek[5][2] =  -Delta/630.;
  Mek[5][3] =  -Delta/420.;
  Mek[5][4] =  2*Delta/315.;
  Mek[5][5] =  4*Delta/315.;
  Mek[5][6] =  2*Delta/315.;
  Mek[5][7] =  Delta/315.;
  Mek[5][8] =  2*Delta/315.;
  Mek[5][9] =  2*Delta/315.;
  Mek[6][0] =  -Delta/630.;
  Mek[6][1] =  -Delta/420.;
  Mek[6][2] =  -Delta/630.;
  Mek[6][3] =  -Delta/420.;
  Mek[6][4] =  2*Delta/315.;
  Mek[6][5] =  2*Delta/315.;
  Mek[6][6] =  4*Delta/315.;
  Mek[6][7] =  2*Delta/315.;
  Mek[6][8] =  Delta/315.;
  Mek[6][9] =  2*Delta/315.;
  Mek[7][0] =  -Delta/630.;
  Mek[7][1] =  -Delta/420.;
  Mek[7][2] =  -Delta/420.;
  Mek[7][3] =  -Delta/630.;
  Mek[7][4] =  2*Delta/315.;
  Mek[7][5] =  Delta/315.;
  Mek[7][6] =  2*Delta/315.;
  Mek[7][7] =  4*Delta/315.;
  Mek[7][8] =  2*Delta/315.;
  Mek[7][9] =  2*Delta/315.;
  Mek[8][0] =  -Delta/420.;
  Mek[8][1] =  -Delta/630.;
  Mek[8][2] =  -Delta/420.;
  Mek[8][3] =  -Delta/630.;
  Mek[8][4] =  2*Delta/315.;
  Mek[8][5] =  2*Delta/315.;
  Mek[8][6] =  Delta/315.;
  Mek[8][7] =  2*Delta/315.;
  Mek[8][8] =  4*Delta/315.;
  Mek[8][9] =  2*Delta/315.;
  Mek[9][0] =  -Delta/420.;
  Mek[9][1] =  -Delta/420.;
  Mek[9][2] =  -Delta/630.;
  Mek[9][3] =  -Delta/630.;
  Mek[9][4] =  Delta/315.;
  Mek[9][5] =  2*Delta/315.;
  Mek[9][6] =  2*Delta/315.;
  Mek[9][7] =  2*Delta/315.;
  Mek[9][8] =  2*Delta/315.;
  Mek[9][9] =  4*Delta/315.;

  return 0;
}


PetscErrorCode Calc3DP1b_Mek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Mek)
{
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscInt n3 = e_inter[en][3];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];
  
  const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 + x1*y0*z2 
            -x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3  
            +x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0  
            +x3*y1*z2 + x3*y2*z0 - x3*y2*z1;

  Mek[0][0] =  7459*Delta/623700.;
  Mek[0][1] =  4523*Delta/1247400.;
  Mek[0][2] =  4523*Delta/1247400.;
  Mek[0][3] =  4523*Delta/1247400.;
  Mek[0][4] =  956*Delta/155925.;
  Mek[1][0] =  4523*Delta/1247400.;
  Mek[1][1] =  7459*Delta/623700.;
  Mek[1][2] =  4523*Delta/1247400.;
  Mek[1][3] =  4523*Delta/1247400.;
  Mek[1][4] =  956*Delta/155925.;
  Mek[2][0] =  4523*Delta/1247400.;
  Mek[2][1] =  4523*Delta/1247400.;
  Mek[2][2] =  7459*Delta/623700.;
  Mek[2][3] =  4523*Delta/1247400.;
  Mek[2][4] =  956*Delta/155925.;
  Mek[3][0] =  4523*Delta/1247400.;
  Mek[3][1] =  4523*Delta/1247400.;
  Mek[3][2] =  4523*Delta/1247400.;
  Mek[3][3] =  7459*Delta/623700.;
  Mek[3][4] =  956*Delta/155925.;
  Mek[4][0] =  956*Delta/155925.;
  Mek[4][1] =  956*Delta/155925.;
  Mek[4][2] =  956*Delta/155925.;
  Mek[4][3] =  956*Delta/155925.;
  Mek[4][4] =  4096*Delta/155925.;

  return 0;
}
#endif