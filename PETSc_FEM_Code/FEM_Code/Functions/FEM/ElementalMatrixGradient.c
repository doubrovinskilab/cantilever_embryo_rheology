#ifndef ElementalMatrixGradient_C
#define ElementalMatrixGradient_C

PetscErrorCode Calc2DP1P1_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];

  Gek_x[0][0] = (y1 - y2)/6.;
  Gek_x[0][1] = (y1 - y2)/6.;
  Gek_x[0][2] = (y1 - y2)/6.;
  Gek_x[1][0] = (y2 - y0)/6.;
  Gek_x[1][1] = (y2 - y0)/6.;
  Gek_x[1][2] = (y2 - y0)/6.;
  Gek_x[2][0] = (y0 - y1)/6.;
  Gek_x[2][1] = (y0 - y1)/6.;
  Gek_x[2][2] = (y0 - y1)/6.;
  
  Gek_y[0][0] = (x2 - x1)/6.;
  Gek_y[0][1] = (x2 - x1)/6.;
  Gek_y[0][2] = (x2 - x1)/6.;
  Gek_y[1][0] = (x0 - x2)/6.;
  Gek_y[1][1] = (x0 - x2)/6.;
  Gek_y[1][2] = (x0 - x2)/6.;
  Gek_y[2][0] = (x1 - x0)/6.;
  Gek_y[2][1] = (x1 - x0)/6.;
  Gek_y[2][2] = (x1 - x0)/6.;

  return 0;
}

PetscErrorCode Calc2DP2P2_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];

  Gek_x[0][0] =  (y1 - y2)/15.;
  Gek_x[0][1] = -(y1 - y2)/30.;
  Gek_x[0][2] = -(y1 - y2)/30.;
  Gek_x[0][3] =  (y1 - y2)/10.;
  Gek_x[0][4] = -(y1 - y2)/30.;
  Gek_x[0][5] =  (y1 - y2)/10.;
  Gek_x[1][0] =  (y0 - y2)/30.;
  Gek_x[1][1] = -(y0 - y2)/15.;
  Gek_x[1][2] =  (y0 - y2)/30.;
  Gek_x[1][3] = -(y0 - y2)/10.;
  Gek_x[1][4] = -(y0 - y2)/10.;
  Gek_x[1][5] =  (y0 - y2)/30.;
  Gek_x[2][0] = -(y0 - y1)/30.;
  Gek_x[2][1] = -(y0 - y1)/30.;
  Gek_x[2][2] =  (y0 - y1)/15.;
  Gek_x[2][3] = -(y0 - y1)/30.;
  Gek_x[2][4] =  (y0 - y1)/10.;
  Gek_x[2][5] =  (y0 - y1)/10.;
  Gek_x[3][0] = -(2.*y0 + y1 - 3.*y2)/30.;
  Gek_x[3][1] =  (y0 + 2.*y1 - 3.*y2)/30.;
  Gek_x[3][2] =  (y0 - y1)/30.;
  Gek_x[3][3] = -4.*(y0 - y1)/15.;
  Gek_x[3][4] = -2.*(y0 - 2.*y1 + y2)/15.;
  Gek_x[3][5] = -2.*(2.*y0 - y1 - y2)/15.;
  Gek_x[4][0] =  (y1 - y2)/30.;
  Gek_x[4][1] =  (3.*y0 - 2.*y1 - y2)/30.;
  Gek_x[4][2] = -(3.*y0 - y1 - 2.*y2)/30.;
  Gek_x[4][3] = 2.*(y0 - 2.*y1 + y2)/15.;
  Gek_x[4][4] = -4.*(y1 - y2)/15.;
  Gek_x[4][5] = -2.*(y0 + y1 - 2.*y2)/15.;
  Gek_x[5][0] =  (2.*y0 - 3.*y1 + y2)/30.;
  Gek_x[5][1] = -(y0 - y2)/30.;
  Gek_x[5][2] = -(y0 - 3.*y1 + 2.*y2)/30.;
  Gek_x[5][3] = 2.*(2.*y0 - y1 - y2)/15.;
  Gek_x[5][4] = 2.*(y0 + y1 - 2.*y2)/15.;
  Gek_x[5][5] = 4.*(y0 - y2)/15.;

  Gek_y[0][0] = -(x1 - x2)/15.;
  Gek_y[0][1] =  (x1 - x2)/30.;
  Gek_y[0][2] =  (x1 - x2)/30.;
  Gek_y[0][3] = -(x1 - x2)/10.;
  Gek_y[0][4] =  (x1 - x2)/30.;
  Gek_y[0][5] = -(x1 - x2)/10.;
  Gek_y[1][0] = -(x0 - x2)/30.;
  Gek_y[1][1] =  (x0 - x2)/15.;
  Gek_y[1][2] = -(x0 - x2)/30.;
  Gek_y[1][3] =  (x0 - x2)/10.;
  Gek_y[1][4] =  (x0 - x2)/10.;
  Gek_y[1][5] = -(x0 - x2)/30.;
  Gek_y[2][0] =  (x0 - x1)/30.;
  Gek_y[2][1] =  (x0 - x1)/30.;
  Gek_y[2][2] = -(x0 - x1)/15.;
  Gek_y[2][3] =  (x0 - x1)/30.;
  Gek_y[2][4] = -(x0 - x1)/10.;
  Gek_y[2][5] = -(x0 - x1)/10.;
  Gek_y[3][0] =  (2.*x0 + x1 - 3*x2)/30.;
  Gek_y[3][1] = -(x0 + 2.*x1 - 3*x2)/30.;
  Gek_y[3][2] = -(x0 - x1)/30.;
  Gek_y[3][3] = 4*(x0 - x1)/15.;
  Gek_y[3][4] = 2.*(x0 - 2.*x1 + x2)/15.;
  Gek_y[3][5] = 2.*(2.*x0 - x1 - x2)/15.;
  Gek_y[4][0] = -(x1 - x2)/30.;
  Gek_y[4][1] = -(3*x0 - 2.*x1 - x2)/30.;
  Gek_y[4][2] =  (3*x0 - x1 - 2.*x2)/30.;
  Gek_y[4][3] = -2.*(x0 - 2.*x1 + x2)/15.;
  Gek_y[4][4] = 4*(x1 - x2)/15.;
  Gek_y[4][5] = 2.*(x0 + x1 - 2.*x2)/15.;
  Gek_y[5][0] = -(2.*x0 - 3*x1 + x2)/30.;
  Gek_y[5][1] =  (x0 - x2)/30.;
  Gek_y[5][2] =  (x0 - 3*x1 + 2.*x2)/30.;
  Gek_y[5][3] = -2.*(2.*x0 - x1 - x2)/15.;
  Gek_y[5][4] = -2.*(x0 + x1 - 2.*x2)/15.;
  Gek_y[5][5] = -4*(x0 - x2)/15.;

	return 0;
}



PetscErrorCode Calc2DP2P1_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];

  Gek_x[0][0] =  (y1 - y2)/6.;
  Gek_x[0][1] =  0.;
  Gek_x[0][2] =  0.;
  Gek_x[1][0] =  0.;
  Gek_x[1][1] = -(y0 - y2)/6.;
  Gek_x[1][2] =  0.;
  Gek_x[2][0] =  0.;
  Gek_x[2][1] =  0.;
  Gek_x[2][2] =  (y0 - y1)/6.;
  Gek_x[3][0] = -(2.*y0 - y1 - y2)/6.;
  Gek_x[3][1] = -(y0 - 2.*y1 + y2)/6.;
  Gek_x[3][2] = -(y0 - y1)/6.;
  Gek_x[4][0] = -(y1 - y2)/6.;
  Gek_x[4][1] =  (y0 - 2.*y1 + y2)/6.;
  Gek_x[4][2] = -(y0 + y1 - 2.*y2)/6.;
  Gek_x[5][0] =  (2.*y0 - y1 - y2)/6.;
  Gek_x[5][1] =  (y0 - y2)/6.;
  Gek_x[5][2] =  (y0 + y1 - 2.*y2)/6.;

  Gek_y[0][0] = -(x1 - x2)/6.;
  Gek_y[0][1] =  0.;
  Gek_y[0][2] =  0.;
  Gek_y[1][0] =  0.;
  Gek_y[1][1] =  (x0 - x2)/6.;
  Gek_y[1][2] =  0.;
  Gek_y[2][0] =  0.;
  Gek_y[2][1] =  0.;
  Gek_y[2][2] = -(x0 - x1)/6.;
  Gek_y[3][0] =  (2.*x0 - x1 - x2)/6.;
  Gek_y[3][1] =  (x0 - 2.*x1 + x2)/6.;
  Gek_y[3][2] =  (x0 - x1)/6.;
  Gek_y[4][0] =  (x1 - x2)/6.;
  Gek_y[4][1] = -(x0 - 2.*x1 + x2)/6.;
  Gek_y[4][2] =  (x0 + x1 - 2.*x2)/6.;
  Gek_y[5][0] = -(2.*x0 - x1 - x2)/6.;
  Gek_y[5][1] = -(x0 - x2)/6.;
  Gek_y[5][2] = -(x0 + x1 - 2.*x2)/6.;

	return 0;
}


PetscErrorCode Calc2DP1P2_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];

  Gek_x[0][0] =  (y1 - y2)/6.;
  Gek_x[0][1] =  0.;
  Gek_x[0][2] =  0.;
  Gek_x[0][3] = -(2.*y0 - y1 - y2)/6.;
  Gek_x[0][4] = -(y1 - y2)/6.;
  Gek_x[0][5] =  (2.*y0 - y1 - y2)/6.;
  Gek_x[1][0] =  0.;
  Gek_x[1][1] = -(y0 - y2)/6.;
  Gek_x[1][2] =  0.;
  Gek_x[1][3] = -(y0 - 2.*y1 + y2)/6.;
  Gek_x[1][4] =  (y0 - 2.*y1 + y2)/6.;
  Gek_x[1][5] =  (y0 - y2)/6.;
  Gek_x[2][0] =  0.;
  Gek_x[2][1] =  0.;
  Gek_x[2][2] =  (y0 - y1)/6.;  
  Gek_x[2][3] = -(y0 - y1)/6.;
  Gek_x[2][4] = -(y0 + y1 - 2.*y2)/6.;
  Gek_x[2][5] =  (y0 + y1 - 2.*y2)/6.;

  Gek_y[0][0] = -(x1 - x2)/6.;
  Gek_y[0][1] =  0.;
  Gek_y[0][2] =  0.;
  Gek_y[0][3] =  (2.*x0 - x1 - x2)/6.;
  Gek_y[0][4] =  (x1 - x2)/6.;
  Gek_y[0][5] = -(2.*x0 - x1 - x2)/6.;
  Gek_y[1][0] =  0.;
  Gek_y[1][1] =  (x0 - x2)/6.;
  Gek_y[1][2] =  0.;
  Gek_y[1][3] =  (x0 - 2.*x1 + x2)/6.;
  Gek_y[1][4] = -(x0 - 2.*x1 + x2)/6.;
  Gek_y[1][5] = -(x0 - x2)/6.;
  Gek_y[2][0] =  0.;
  Gek_y[2][1] =  0.;
  Gek_y[2][2] = -(x0 - x1)/6.;
  Gek_y[2][3] =  (x0 - x1)/6.;
  Gek_y[2][4] =  (x0 + x1 - 2.*x2)/6.;
  Gek_y[2][5] = -(x0 + x1 - 2.*x2)/6.;

	return 0;
}


PetscErrorCode Calc2DP1bP1_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];

  Gek_x[0][0] =  29.*(y1 - y2)/120.;
  Gek_x[0][1] =  -(9.*y0 - 20.*y1 + 11.*y2)/120.;
  Gek_x[0][2] =   (9.*y0 + 11.*y1 - 20.*y2)/120.;
  Gek_x[1][0] =  -(20.*y0 - 9.*y1 - 11.*y2)/120.;
  Gek_x[1][1] =  -29.*(y0 - y2)/120.;
  Gek_x[1][2] =  -(11.*y0 + 9.*y1 - 20.*y2)/120.;
  Gek_x[2][0] =   (20.*y0 - 11.*y1 - 9.*y2)/120.;
  Gek_x[2][1] =   (11.*y0 - 20.*y1 + 9.*y2)/120.;
  Gek_x[2][2] =  29.*(y0 - y1)/120.;
  Gek_x[3][0] =  -9.*(y1 - y2)/40.;
  Gek_x[3][1] =   9.*(y0 - y2)/40.;
  Gek_x[3][2] =  -9.*(y0 - y1)/40.;
  
  Gek_y[0][0] =  -29.*(x1 - x2)/120.;
  Gek_y[0][1] =   (9.*x0 - 20.*x1 + 11.*x2)/120.;
  Gek_y[0][2] =  -(9.*x0 + 11.*x1 - 20.*x2)/120.;
  Gek_y[1][0] =   (20.*x0 - 9.*x1 - 11.*x2)/120.;
  Gek_y[1][1] =  29.*(x0 - x2)/120.;
  Gek_y[1][2] =   (11.*x0 + 9.*x1 - 20.*x2)/120.;
  Gek_y[2][0] =  -(20.*x0 - 11.*x1 - 9.*x2)/120.;
  Gek_y[2][1] =  -(11.*x0 - 20.*x1 + 9.*x2)/120.;
  Gek_y[2][2] = -29.*(x0 - x1)/120.;
  Gek_y[3][0] =   9.*(x1 - x2)/40.;
  Gek_y[3][1] =  -9.*(x0 - x2)/40.;
  Gek_y[3][2] =   9.*(x0 - x1)/40.;

  return 0;
}

PetscErrorCode Calc2DP1P1b_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{ 
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];

  Gek_x[0][0] =  29.*(y1 - y2)/120.;
  Gek_x[0][1] =  -(20.*y0 - 9.*y1 - 11.*y2)/120.;
  Gek_x[0][2] =   (20.*y0 - 11.*y1 - 9.*y2)/120.;
  Gek_x[0][3] =  -9.*(y1 - y2)/40.;
  Gek_x[1][0] =  -(9.*y0 - 20.*y1 + 11.*y2)/120.;
  Gek_x[1][1] =  -29.*(y0 - y2)/120.;
  Gek_x[1][2] =   (11.*y0 - 20.*y1 + 9.*y2)/120.;
  Gek_x[1][3] =   9.*(y0 - y2)/40.;
  Gek_x[2][0] =   (9.*y0 + 11.*y1 - 20.*y2)/120.;
  Gek_x[2][1] =  -(11.*y0 + 9.*y1 - 20.*y2)/120.;
  Gek_x[2][2] =  29.*(y0 - y1)/120.;
  Gek_x[2][3] =  -9.*(y0 - y1)/40.;
  
  Gek_y[0][0] =  -29.*(x1 - x2)/120.;
  Gek_y[0][1] =   (20.*x0 - 9.*x1 - 11.*x2)/120.;
  Gek_y[0][2] =  -(20.*x0 - 11.*x1 - 9.*x2)/120.;
  Gek_y[0][3] =   9.*(x1 - x2)/40.;
  Gek_y[1][0] =   (9.*x0 - 20.*x1 + 11.*x2)/120.;
  Gek_y[1][1] =  29.*(x0 - x2)/120.;
  Gek_y[1][2] =  -(11.*x0 - 20.*x1 + 9.*x2)/120.;
  Gek_y[1][3] =  -9.*(x0 - x2)/40.;
  Gek_y[2][0] =  -(9.*x0 + 11.*x1 - 20.*x2)/120.;
  Gek_y[2][1] =   (11.*x0 + 9.*x1 - 20.*x2)/120.;
  Gek_y[2][2] = -29.*(x0 - x1)/120.;
  Gek_y[2][3] =   9.*(x0 - x1)/40.;

  return 0;
}


PetscErrorCode Calc2DP1bP1b_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];

  Gek_x[0][0] =  (y1 - y2)/6.;
  Gek_x[0][1] = -(9.*y0 - 11.*y1 +  2.*y2)/120.;
  Gek_x[0][2] =  (9.*y0 +  2.*y1 - 11.*y2)/120.;
  Gek_x[0][3] =   9.*(y1 - y2)/40.;
  Gek_x[1][0] = -(11.*y0 - 9.*y1 - 2.*y2)/120.;
  Gek_x[1][1] = -(y0 - y2)/6.;
  Gek_x[1][2] = -(2.*y0 + 9.*y1 - 11.*y2)/120.;
  Gek_x[1][3] = -9.*(y0 - y2)/40.;
  Gek_x[2][0] =  (11.*y0 - 2.*y1 - 9.*y2)/120.;
  Gek_x[2][1] =  (2.*y0 - 11.*y1 + 9.*y2)/120.;
  Gek_x[2][2] =     (y0 - y1)/6.;
  Gek_x[2][3] =  9.*(y0 - y1)/40.;
  Gek_x[3][0] = -9.*(y1 - y2)/40.;
  Gek_x[3][1] =  9.*(y0 - y2)/40.;
  Gek_x[3][2] = -9.*(y0 - y1)/40.;
  Gek_x[3][3] =  0.;
  
  Gek_y[0][0] =  -(x1 - x2)/6.;
  Gek_y[0][1] =   (9.*x0 - 11.*x1 + 2.*x2)/120.;
  Gek_y[0][2] =  -(9.*x0 + 2.*x1 - 11.*x2)/120.;
  Gek_y[0][3] =  -9.*(x1 - x2)/40.;
  Gek_y[1][0] =  (11.*x0 - 9.*x1 - 2.*x2)/120.;
  Gek_y[1][1] =  (x0 - x2)/6.;
  Gek_y[1][2] =  (2.*x0 + 9.*x1 - 11.*x2)/120.;
  Gek_y[1][3] =   9.*(x0 - x2)/40.;
  Gek_y[2][0] =  -(11.*x0 - 2.*x1 - 9.*x2)/120.;
  Gek_y[2][1] =  -(2.*x0 - 11.*x1 + 9.*x2)/120.;
  Gek_y[2][2] =  -(x0 - x1)/6.;
  Gek_y[2][3] =  -9.*(x0 - x1)/40.;
  Gek_y[3][0] =   9.*(x1 - x2)/40.;
  Gek_y[3][1] =  -9.*(x0 - x2)/40.;
  Gek_y[3][2] =   9.*(x0 - x1)/40.;
  Gek_y[3][3] =   0.;

  return 0;
}


PetscErrorCode Calc3DP1P1_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscInt n3 = e_inter[en][3];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];

  Gek_x[0][0] = -(y1*z2 - y1*z3 - y2*z1 + y2*z3 + y3*z1 - y3*z2)/4.0;
  Gek_x[0][1] = -(y1*z2 - y1*z3 - y2*z1 + y2*z3 + y3*z1 - y3*z2)/4.0;
  Gek_x[0][2] = -(y1*z2 - y1*z3 - y2*z1 + y2*z3 + y3*z1 - y3*z2)/4.0;
  Gek_x[0][3] = -(y1*z2 - y1*z3 - y2*z1 + y2*z3 + y3*z1 - y3*z2)/4.0;
  Gek_x[1][0] =  (y0*z2 - y0*z3 - y2*z0 + y2*z3 + y3*z0 - y3*z2)/4.0;
  Gek_x[1][1] =  (y0*z2 - y0*z3 - y2*z0 + y2*z3 + y3*z0 - y3*z2)/4.0;
  Gek_x[1][2] =  (y0*z2 - y0*z3 - y2*z0 + y2*z3 + y3*z0 - y3*z2)/4.0;
  Gek_x[1][3] =  (y0*z2 - y0*z3 - y2*z0 + y2*z3 + y3*z0 - y3*z2)/4.0;
  Gek_x[2][0] = -(y0*z1 - y0*z3 - y1*z0 + y1*z3 + y3*z0 - y3*z1)/4.0;
  Gek_x[2][1] = -(y0*z1 - y0*z3 - y1*z0 + y1*z3 + y3*z0 - y3*z1)/4.0;
  Gek_x[2][2] = -(y0*z1 - y0*z3 - y1*z0 + y1*z3 + y3*z0 - y3*z1)/4.0;
  Gek_x[2][3] = -(y0*z1 - y0*z3 - y1*z0 + y1*z3 + y3*z0 - y3*z1)/4.0;
  Gek_x[3][0] =  (y0*z1 - y0*z2 - y1*z0 + y1*z2 + y2*z0 - y2*z1)/4.0;
  Gek_x[3][1] =  (y0*z1 - y0*z2 - y1*z0 + y1*z2 + y2*z0 - y2*z1)/4.0;
  Gek_x[3][2] =  (y0*z1 - y0*z2 - y1*z0 + y1*z2 + y2*z0 - y2*z1)/4.0;
  Gek_x[3][3] =  (y0*z1 - y0*z2 - y1*z0 + y1*z2 + y2*z0 - y2*z1)/4.0;
  
  Gek_y[0][0] =  (x1*z2 - x1*z3 - x2*z1 + x2*z3 + x3*z1 - x3*z2)/4.0;
  Gek_y[0][1] =  (x1*z2 - x1*z3 - x2*z1 + x2*z3 + x3*z1 - x3*z2)/4.0;
  Gek_y[0][2] =  (x1*z2 - x1*z3 - x2*z1 + x2*z3 + x3*z1 - x3*z2)/4.0;
  Gek_y[0][3] =  (x1*z2 - x1*z3 - x2*z1 + x2*z3 + x3*z1 - x3*z2)/4.0;
  Gek_y[1][0] = -(x0*z2 - x0*z3 - x2*z0 + x2*z3 + x3*z0 - x3*z2)/4.0;
  Gek_y[1][1] = -(x0*z2 - x0*z3 - x2*z0 + x2*z3 + x3*z0 - x3*z2)/4.0;
  Gek_y[1][2] = -(x0*z2 - x0*z3 - x2*z0 + x2*z3 + x3*z0 - x3*z2)/4.0;
  Gek_y[1][3] = -(x0*z2 - x0*z3 - x2*z0 + x2*z3 + x3*z0 - x3*z2)/4.0;
  Gek_y[2][0] =  (x0*z1 - x0*z3 - x1*z0 + x1*z3 + x3*z0 - x3*z1)/4.0;
  Gek_y[2][1] =  (x0*z1 - x0*z3 - x1*z0 + x1*z3 + x3*z0 - x3*z1)/4.0;
  Gek_y[2][2] =  (x0*z1 - x0*z3 - x1*z0 + x1*z3 + x3*z0 - x3*z1)/4.0;
  Gek_y[2][3] =  (x0*z1 - x0*z3 - x1*z0 + x1*z3 + x3*z0 - x3*z1)/4.0;
  Gek_y[3][0] = -(x0*z1 - x0*z2 - x1*z0 + x1*z2 + x2*z0 - x2*z1)/4.0;
  Gek_y[3][1] = -(x0*z1 - x0*z2 - x1*z0 + x1*z2 + x2*z0 - x2*z1)/4.0;
  Gek_y[3][2] = -(x0*z1 - x0*z2 - x1*z0 + x1*z2 + x2*z0 - x2*z1)/4.0;
  Gek_y[3][3] = -(x0*z1 - x0*z2 - x1*z0 + x1*z2 + x2*z0 - x2*z1)/4.0;
  
  Gek_z[0][0] = -(x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2)/4.0;
  Gek_z[0][1] = -(x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2)/4.0;
  Gek_z[0][2] = -(x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2)/4.0;
  Gek_z[0][3] = -(x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2)/4.0;
  Gek_z[1][0] =  (x0*y2 - x0*y3 - x2*y0 + x2*y3 + x3*y0 - x3*y2)/4.0;
  Gek_z[1][1] =  (x0*y2 - x0*y3 - x2*y0 + x2*y3 + x3*y0 - x3*y2)/4.0;
  Gek_z[1][2] =  (x0*y2 - x0*y3 - x2*y0 + x2*y3 + x3*y0 - x3*y2)/4.0;
  Gek_z[1][3] =  (x0*y2 - x0*y3 - x2*y0 + x2*y3 + x3*y0 - x3*y2)/4.0;
  Gek_z[2][0] = -(x0*y1 - x0*y3 - x1*y0 + x1*y3 + x3*y0 - x3*y1)/4.0;
  Gek_z[2][1] = -(x0*y1 - x0*y3 - x1*y0 + x1*y3 + x3*y0 - x3*y1)/4.0;
  Gek_z[2][2] = -(x0*y1 - x0*y3 - x1*y0 + x1*y3 + x3*y0 - x3*y1)/4.0;
  Gek_z[2][3] = -(x0*y1 - x0*y3 - x1*y0 + x1*y3 + x3*y0 - x3*y1)/4.0;
  Gek_z[3][0] =  (x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)/4.0;
  Gek_z[3][1] =  (x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)/4.0;
  Gek_z[3][2] =  (x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)/4.0;
  Gek_z[3][3] =  (x0*y1 - x0*y2 - x1*y0 + x1*y2 + x2*y0 - x2*y1)/4.0;

  return 0;
}

PetscErrorCode Calc3DP2P2_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscInt n3 = e_inter[en][3];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];

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
  
  Gek_x[0][0] = -C0yz/120.;
  Gek_x[0][1] = C0yz/360.;
  Gek_x[0][2] = C0yz/360.;
  Gek_x[0][3] = C0yz/360.;
  Gek_x[0][4] = -C0yz/90.;
  Gek_x[0][5] = C0yz/90.;
  Gek_x[0][6] = -C0yz/90.;
  Gek_x[0][7] = -C0yz/90.;
  Gek_x[0][8] = C0yz/90.;
  Gek_x[0][9] = C0yz/90.;
  Gek_x[1][0] = C1yz/360.;
  Gek_x[1][1] = -C1yz/120.;
  Gek_x[1][2] = C1yz/360.;
  Gek_x[1][3] = C1yz/360.;
  Gek_x[1][4] = -C1yz/90.;
  Gek_x[1][5] = -C1yz/90.;
  Gek_x[1][6] = C1yz/90.;
  Gek_x[1][7] = C1yz/90.;
  Gek_x[1][8] = -C1yz/90.;
  Gek_x[1][9] = C1yz/90.;
  Gek_x[2][0] = C2yz/360.;
  Gek_x[2][1] = C2yz/360.;
  Gek_x[2][2] = -C2yz/120.;
  Gek_x[2][3] = C2yz/360.;
  Gek_x[2][4] = C2yz/90.;
  Gek_x[2][5] = -C2yz/90.;
  Gek_x[2][6] = -C2yz/90.;
  Gek_x[2][7] = C2yz/90.;
  Gek_x[2][8] = C2yz/90.;
  Gek_x[2][9] = -C2yz/90.;
  Gek_x[3][0] = C3yz/360.;
  Gek_x[3][1] = C3yz/360.;
  Gek_x[3][2] = C3yz/360.;
  Gek_x[3][3] = -C3yz/120.;
  Gek_x[3][4] = C3yz/90.;
  Gek_x[3][5] = C3yz/90.;
  Gek_x[3][6] = C3yz/90.;
  Gek_x[3][7] = -C3yz/90.;
  Gek_x[3][8] = -C3yz/90.;
  Gek_x[3][9] = -C3yz/90.;
  Gek_x[4][0] = C0yz/90.;
  Gek_x[4][1] = C1yz/90.;
  Gek_x[4][2] = (C0yz + C1yz)/90.;
  Gek_x[4][3] = (C0yz + C1yz)/90.;
  Gek_x[4][4] = -2.*(C0yz + C1yz)/45.;
  Gek_x[4][5] = -(2.*C0yz + C1yz)/45.;
  Gek_x[4][6] = -(C0yz + 2.*C1yz)/45.;
  Gek_x[4][7] = -(C0yz + 2.*C1yz)/45.;
  Gek_x[4][8] = -(2.*C0yz + C1yz)/45.;
  Gek_x[4][9] = -(C0yz + C1yz)/45.;
  Gek_x[5][0] = (C1yz + C2yz)/90.;
  Gek_x[5][1] = C1yz/90.;
  Gek_x[5][2] = C2yz/90.;
  Gek_x[5][3] = (C1yz + C2yz)/90.;
  Gek_x[5][4] = -(C1yz + 2.*C2yz)/45.;
  Gek_x[5][5] = -2.*(C1yz + C2yz)/45.;
  Gek_x[5][6] = -(2.*C1yz + C2yz)/45.;
  Gek_x[5][7] = -(C1yz + C2yz)/45.;
  Gek_x[5][8] = -(C1yz + 2.*C2yz)/45.;
  Gek_x[5][9] = -(2.*C1yz + C2yz)/45.;
  Gek_x[6][0] = C0yz/90.;
  Gek_x[6][1] = (C0yz + C2yz)/90.;
  Gek_x[6][2] = C2yz/90.;
  Gek_x[6][3] = (C0yz + C2yz)/90.;
  Gek_x[6][4] = -(C0yz + 2.*C2yz)/45.;
  Gek_x[6][5] = -(2.*C0yz + C2yz)/45.;
  Gek_x[6][6] = -2.*(C0yz + C2yz)/45.;
  Gek_x[6][7] = -(C0yz + 2.*C2yz)/45.;
  Gek_x[6][8] = -(C0yz + C2yz)/45.;
  Gek_x[6][9] = -(2.*C0yz + C2yz)/45.;
  Gek_x[7][0] = C0yz/90.;
  Gek_x[7][1] = (C0yz + C3yz)/90.;
  Gek_x[7][2] = (C0yz + C3yz)/90.;
  Gek_x[7][3] = C3yz/90.;
  Gek_x[7][4] = -(C0yz + 2.*C3yz)/45.;
  Gek_x[7][5] = -(C0yz + C3yz)/45.;
  Gek_x[7][6] = -(C0yz + 2.*C3yz)/45.;
  Gek_x[7][7] = -2.*(C0yz + C3yz)/45.;
  Gek_x[7][8] = -(2.*C0yz + C3yz)/45.;
  Gek_x[7][9] = -(2.*C0yz + C3yz)/45.;
  Gek_x[8][0] = (C1yz + C3yz)/90.;
  Gek_x[8][1] = C1yz/90.;
  Gek_x[8][2] = (C1yz + C3yz)/90.;
  Gek_x[8][3] = C3yz/90.;
  Gek_x[8][4] = -(C1yz + 2.*C3yz)/45.;
  Gek_x[8][5] = -(C1yz + 2.*C3yz)/45.;
  Gek_x[8][6] = -(C1yz + C3yz)/45.;
  Gek_x[8][7] = -(2.*C1yz + C3yz)/45.;
  Gek_x[8][8] = -2.*(C1yz + C3yz)/45.;
  Gek_x[8][9] = -(2.*C1yz + C3yz)/45.;
  Gek_x[9][0] = (C2yz + C3yz)/90.;
  Gek_x[9][1] = (C2yz + C3yz)/90.;
  Gek_x[9][2] = C2yz/90.;
  Gek_x[9][3] = C3yz/90.;
  Gek_x[9][4] = -(C2yz + C3yz)/45.;
  Gek_x[9][5] = -(C2yz + 2.*C3yz)/45.;
  Gek_x[9][6] = -(C2yz + 2.*C3yz)/45.;
  Gek_x[9][7] = -(2.*C2yz + C3yz)/45.;
  Gek_x[9][8] = -(2.*C2yz + C3yz)/45.;
  Gek_x[9][9] = -2.*(C2yz + C3yz)/45.;

  Gek_y[0][0] = -C0xz/120.;
  Gek_y[0][1] = C0xz/360.;
  Gek_y[0][2] = C0xz/360.;
  Gek_y[0][3] = C0xz/360.;
  Gek_y[0][4] = -C0xz/90.;
  Gek_y[0][5] = C0xz/90.;
  Gek_y[0][6] = -C0xz/90.;
  Gek_y[0][7] = -C0xz/90.;
  Gek_y[0][8] = C0xz/90.;
  Gek_y[0][9] = C0xz/90.;
  Gek_y[1][0] = C1xz/360.;
  Gek_y[1][1] = -C1xz/120.;
  Gek_y[1][2] = C1xz/360.;
  Gek_y[1][3] = C1xz/360.;
  Gek_y[1][4] = -C1xz/90.;
  Gek_y[1][5] = -C1xz/90.;
  Gek_y[1][6] = C1xz/90.;
  Gek_y[1][7] = C1xz/90.;
  Gek_y[1][8] = -C1xz/90.;
  Gek_y[1][9] = C1xz/90.;
  Gek_y[2][0] = C2xz/360.;
  Gek_y[2][1] = C2xz/360.;
  Gek_y[2][2] = -C2xz/120.;
  Gek_y[2][3] = C2xz/360.;
  Gek_y[2][4] = C2xz/90.;
  Gek_y[2][5] = -C2xz/90.;
  Gek_y[2][6] = -C2xz/90.;
  Gek_y[2][7] = C2xz/90.;
  Gek_y[2][8] = C2xz/90.;
  Gek_y[2][9] = -C2xz/90.;
  Gek_y[3][0] = C3xz/360.;
  Gek_y[3][1] = C3xz/360.;
  Gek_y[3][2] = C3xz/360.;
  Gek_y[3][3] = -C3xz/120.;
  Gek_y[3][4] = C3xz/90.;
  Gek_y[3][5] = C3xz/90.;
  Gek_y[3][6] = C3xz/90.;
  Gek_y[3][7] = -C3xz/90.;
  Gek_y[3][8] = -C3xz/90.;
  Gek_y[3][9] = -C3xz/90.;
  Gek_y[4][0] = C0xz/90.;
  Gek_y[4][1] = C1xz/90.;
  Gek_y[4][2] = (C0xz + C1xz)/90.;
  Gek_y[4][3] = (C0xz + C1xz)/90.;
  Gek_y[4][4] = -2.*(C0xz + C1xz)/45.;
  Gek_y[4][5] = -(2.*C0xz + C1xz)/45.;
  Gek_y[4][6] = -(C0xz + 2.*C1xz)/45.;
  Gek_y[4][7] = -(C0xz + 2.*C1xz)/45.;
  Gek_y[4][8] = -(2.*C0xz + C1xz)/45.;
  Gek_y[4][9] = -(C0xz + C1xz)/45.;
  Gek_y[5][0] = (C1xz + C2xz)/90.;
  Gek_y[5][1] = C1xz/90.;
  Gek_y[5][2] = C2xz/90.;
  Gek_y[5][3] = (C1xz + C2xz)/90.;
  Gek_y[5][4] = -(C1xz + 2.*C2xz)/45.;
  Gek_y[5][5] = -2.*(C1xz + C2xz)/45.;
  Gek_y[5][6] = -(2.*C1xz + C2xz)/45.;
  Gek_y[5][7] = -(C1xz + C2xz)/45.;
  Gek_y[5][8] = -(C1xz + 2.*C2xz)/45.;
  Gek_y[5][9] = -(2.*C1xz + C2xz)/45.;
  Gek_y[6][0] = C0xz/90.;
  Gek_y[6][1] = (C0xz + C2xz)/90.;
  Gek_y[6][2] = C2xz/90.;
  Gek_y[6][3] = (C0xz + C2xz)/90.;
  Gek_y[6][4] = -(C0xz + 2.*C2xz)/45.;
  Gek_y[6][5] = -(2.*C0xz + C2xz)/45.;
  Gek_y[6][6] = -2.*(C0xz + C2xz)/45.;
  Gek_y[6][7] = -(C0xz + 2.*C2xz)/45.;
  Gek_y[6][8] = -(C0xz + C2xz)/45.;
  Gek_y[6][9] = -(2.*C0xz + C2xz)/45.;
  Gek_y[7][0] = C0xz/90.;
  Gek_y[7][1] = (C0xz + C3xz)/90.;
  Gek_y[7][2] = (C0xz + C3xz)/90.;
  Gek_y[7][3] = C3xz/90.;
  Gek_y[7][4] = -(C0xz + 2.*C3xz)/45.;
  Gek_y[7][5] = -(C0xz + C3xz)/45.;
  Gek_y[7][6] = -(C0xz + 2.*C3xz)/45.;
  Gek_y[7][7] = -2.*(C0xz + C3xz)/45.;
  Gek_y[7][8] = -(2.*C0xz + C3xz)/45.;
  Gek_y[7][9] = -(2.*C0xz + C3xz)/45.;
  Gek_y[8][0] = (C1xz + C3xz)/90.;
  Gek_y[8][1] = C1xz/90.;
  Gek_y[8][2] = (C1xz + C3xz)/90.;
  Gek_y[8][3] = C3xz/90.;
  Gek_y[8][4] = -(C1xz + 2.*C3xz)/45.;
  Gek_y[8][5] = -(C1xz + 2.*C3xz)/45.;
  Gek_y[8][6] = -(C1xz + C3xz)/45.;
  Gek_y[8][7] = -(2.*C1xz + C3xz)/45.;
  Gek_y[8][8] = -2.*(C1xz + C3xz)/45.;
  Gek_y[8][9] = -(2.*C1xz + C3xz)/45.;
  Gek_y[9][0] = (C2xz + C3xz)/90.;
  Gek_y[9][1] = (C2xz + C3xz)/90.;
  Gek_y[9][2] = C2xz/90.;
  Gek_y[9][3] = C3xz/90.;
  Gek_y[9][4] = -(C2xz + C3xz)/45.;
  Gek_y[9][5] = -(C2xz + 2.*C3xz)/45.;
  Gek_y[9][6] = -(C2xz + 2.*C3xz)/45.;
  Gek_y[9][7] = -(2.*C2xz + C3xz)/45.;
  Gek_y[9][8] = -(2.*C2xz + C3xz)/45.;
  Gek_y[9][9] = -2.*(C2xz + C3xz)/45.;
  
  Gek_z[0][0] = -C0xy/120.;
  Gek_z[0][1] = C0xy/360.;
  Gek_z[0][2] = C0xy/360.;
  Gek_z[0][3] = C0xy/360.;
  Gek_z[0][4] = -C0xy/90.;
  Gek_z[0][5] = C0xy/90.;
  Gek_z[0][6] = -C0xy/90.;
  Gek_z[0][7] = -C0xy/90.;
  Gek_z[0][8] = C0xy/90.;
  Gek_z[0][9] = C0xy/90.;
  Gek_z[1][0] = C1xy/360.;
  Gek_z[1][1] = -C1xy/120.;
  Gek_z[1][2] = C1xy/360.;
  Gek_z[1][3] = C1xy/360.;
  Gek_z[1][4] = -C1xy/90.;
  Gek_z[1][5] = -C1xy/90.;
  Gek_z[1][6] = C1xy/90.;
  Gek_z[1][7] = C1xy/90.;
  Gek_z[1][8] = -C1xy/90.;
  Gek_z[1][9] = C1xy/90.;
  Gek_z[2][0] = C2xy/360.;
  Gek_z[2][1] = C2xy/360.;
  Gek_z[2][2] = -C2xy/120.;
  Gek_z[2][3] = C2xy/360.;
  Gek_z[2][4] = C2xy/90.;
  Gek_z[2][5] = -C2xy/90.;
  Gek_z[2][6] = -C2xy/90.;
  Gek_z[2][7] = C2xy/90.;
  Gek_z[2][8] = C2xy/90.;
  Gek_z[2][9] = -C2xy/90.;
  Gek_z[3][0] = C3xy/360.;
  Gek_z[3][1] = C3xy/360.;
  Gek_z[3][2] = C3xy/360.;
  Gek_z[3][3] = -C3xy/120.;
  Gek_z[3][4] = C3xy/90.;
  Gek_z[3][5] = C3xy/90.;
  Gek_z[3][6] = C3xy/90.;
  Gek_z[3][7] = -C3xy/90.;
  Gek_z[3][8] = -C3xy/90.;
  Gek_z[3][9] = -C3xy/90.;
  Gek_z[4][0] = C0xy/90.;
  Gek_z[4][1] = C1xy/90.;
  Gek_z[4][2] = (C0xy + C1xy)/90.;
  Gek_z[4][3] = (C0xy + C1xy)/90.;
  Gek_z[4][4] = -2.*(C0xy + C1xy)/45.;
  Gek_z[4][5] = -(2.*C0xy + C1xy)/45.;
  Gek_z[4][6] = -(C0xy + 2.*C1xy)/45.;
  Gek_z[4][7] = -(C0xy + 2.*C1xy)/45.;
  Gek_z[4][8] = -(2.*C0xy + C1xy)/45.;
  Gek_z[4][9] = -(C0xy + C1xy)/45.;
  Gek_z[5][0] = (C1xy + C2xy)/90.;
  Gek_z[5][1] = C1xy/90.;
  Gek_z[5][2] = C2xy/90.;
  Gek_z[5][3] = (C1xy + C2xy)/90.;
  Gek_z[5][4] = -(C1xy + 2.*C2xy)/45.;
  Gek_z[5][5] = -2.*(C1xy + C2xy)/45.;
  Gek_z[5][6] = -(2.*C1xy + C2xy)/45.;
  Gek_z[5][7] = -(C1xy + C2xy)/45.;
  Gek_z[5][8] = -(C1xy + 2.*C2xy)/45.;
  Gek_z[5][9] = -(2.*C1xy + C2xy)/45.;
  Gek_z[6][0] = C0xy/90.;
  Gek_z[6][1] = (C0xy + C2xy)/90.;
  Gek_z[6][2] = C2xy/90.;
  Gek_z[6][3] = (C0xy + C2xy)/90.;
  Gek_z[6][4] = -(C0xy + 2.*C2xy)/45.;
  Gek_z[6][5] = -(2.*C0xy + C2xy)/45.;
  Gek_z[6][6] = -2.*(C0xy + C2xy)/45.;
  Gek_z[6][7] = -(C0xy + 2.*C2xy)/45.;
  Gek_z[6][8] = -(C0xy + C2xy)/45.;
  Gek_z[6][9] = -(2.*C0xy + C2xy)/45.;
  Gek_z[7][0] = C0xy/90.;
  Gek_z[7][1] = (C0xy + C3xy)/90.;
  Gek_z[7][2] = (C0xy + C3xy)/90.;
  Gek_z[7][3] = C3xy/90.;
  Gek_z[7][4] = -(C0xy + 2.*C3xy)/45.;
  Gek_z[7][5] = -(C0xy + C3xy)/45.;
  Gek_z[7][6] = -(C0xy + 2.*C3xy)/45.;
  Gek_z[7][7] = -2.*(C0xy + C3xy)/45.;
  Gek_z[7][8] = -(2.*C0xy + C3xy)/45.;
  Gek_z[7][9] = -(2.*C0xy + C3xy)/45.;
  Gek_z[8][0] = (C1xy + C3xy)/90.;
  Gek_z[8][1] = C1xy/90.;
  Gek_z[8][2] = (C1xy + C3xy)/90.;
  Gek_z[8][3] = C3xy/90.;
  Gek_z[8][4] = -(C1xy + 2.*C3xy)/45.;
  Gek_z[8][5] = -(C1xy + 2.*C3xy)/45.;
  Gek_z[8][6] = -(C1xy + C3xy)/45.;
  Gek_z[8][7] = -(2.*C1xy + C3xy)/45.;
  Gek_z[8][8] = -2.*(C1xy + C3xy)/45.;
  Gek_z[8][9] = -(2.*C1xy + C3xy)/45.;
  Gek_z[9][0] = (C2xy + C3xy)/90.;
  Gek_z[9][1] = (C2xy + C3xy)/90.;
  Gek_z[9][2] = C2xy/90.;
  Gek_z[9][3] = C3xy/90.;
  Gek_z[9][4] = -(C2xy + C3xy)/45.;
  Gek_z[9][5] = -(C2xy + 2.*C3xy)/45.;
  Gek_z[9][6] = -(C2xy + 2.*C3xy)/45.;
  Gek_z[9][7] = -(2.*C2xy + C3xy)/45.;
  Gek_z[9][8] = -(2.*C2xy + C3xy)/45.;
  Gek_z[9][9] = -2.*(C2xy + C3xy)/45.;

	return 0;
}


PetscErrorCode Calc3DP2P1_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscInt n3 = e_inter[en][3];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];

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

  Gek_x[0][0] = -C0yz/40.;
  Gek_x[0][1] =  C0yz/120.;
  Gek_x[0][2] =  C0yz/120.;
  Gek_x[0][3] =  C0yz/120.;
  Gek_x[1][0] =  C1yz/120.;
  Gek_x[1][1] = -C1yz/40.;
  Gek_x[1][2] =  C1yz/120.;
  Gek_x[1][3] =  C1yz/120.;
  Gek_x[2][0] =  C2yz/120.;
  Gek_x[2][1] =  C2yz/120.;
  Gek_x[2][2] = -C2yz/40.;
  Gek_x[2][3] =  C2yz/120.;
  Gek_x[3][0] =  C3yz/120.;
  Gek_x[3][1] =  C3yz/120.;
  Gek_x[3][2] =  C3yz/120.;
  Gek_x[3][3] = -C3yz/40.;
  Gek_x[4][0] = -(C0yz + 2.*C1yz)/30.;
  Gek_x[4][1] = -(2.*C0yz + C1yz)/30.;
  Gek_x[4][2] = -(C0yz + C1yz)/30.;
  Gek_x[4][3] = -(C0yz + C1yz)/30.;
  Gek_x[5][0] = -(C1yz + C2yz)/30.;
  Gek_x[5][1] = -(C1yz + 2.*C2yz)/30.;
  Gek_x[5][2] = -(2.*C1yz + C2yz)/30.;
  Gek_x[5][3] = -(C1yz + C2yz)/30.;
  Gek_x[6][0] = -(C0yz + 2.*C2yz)/30.;
  Gek_x[6][1] = -(C0yz + C2yz)/30.;
  Gek_x[6][2] = -(2.*C0yz + C2yz)/30.;
  Gek_x[6][3] = -(C0yz + C2yz)/30.;
  Gek_x[7][0] = -(C0yz + 2.*C3yz)/30.;
  Gek_x[7][1] = -(C0yz + C3yz)/30.;
  Gek_x[7][2] = -(C0yz + C3yz)/30.;
  Gek_x[7][3] = -(2.*C0yz + C3yz)/30.;
  Gek_x[8][0] = -(C1yz + C3yz)/30.;
  Gek_x[8][1] = -(C1yz + 2.*C3yz)/30.;
  Gek_x[8][2] = -(C1yz + C3yz)/30.;
  Gek_x[8][3] = -(2.*C1yz + C3yz)/30.;
  Gek_x[9][0] = -(C2yz + C3yz)/30.;
  Gek_x[9][1] = -(C2yz + C3yz)/30.;
  Gek_x[9][2] = -(C2yz + 2.*C3yz)/30.;
  Gek_x[9][3] = -(2.*C2yz + C3yz)/30.;
  
  Gek_y[0][0] = -C0xz/40.;
  Gek_y[0][1] = C0xz/120.;
  Gek_y[0][2] = C0xz/120.;
  Gek_y[0][3] = C0xz/120.;
  Gek_y[1][0] = C1xz/120.;
  Gek_y[1][1] = -C1xz/40.;
  Gek_y[1][2] = C1xz/120.;
  Gek_y[1][3] = C1xz/120.;
  Gek_y[2][0] = C2xz/120.;
  Gek_y[2][1] = C2xz/120.;
  Gek_y[2][2] = -C2xz/40.;
  Gek_y[2][3] = C2xz/120.;
  Gek_y[3][0] = C3xz/120.;
  Gek_y[3][1] = C3xz/120.;
  Gek_y[3][2] = C3xz/120.;
  Gek_y[3][3] = -C3xz/40.;
  Gek_y[4][0] = -(C0xz + 2.*C1xz)/30.;
  Gek_y[4][1] = -(2.*C0xz + C1xz)/30.;
  Gek_y[4][2] = -(C0xz + C1xz)/30.;
  Gek_y[4][3] = -(C0xz + C1xz)/30.;
  Gek_y[5][0] = -(C1xz + C2xz)/30.;
  Gek_y[5][1] = -(C1xz + 2.*C2xz)/30.;
  Gek_y[5][2] = -(2.*C1xz + C2xz)/30.;
  Gek_y[5][3] = -(C1xz + C2xz)/30.;
  Gek_y[6][0] = -(C0xz + 2.*C2xz)/30.;
  Gek_y[6][1] = -(C0xz + C2xz)/30.;
  Gek_y[6][2] = -(2.*C0xz + C2xz)/30.;
  Gek_y[6][3] = -(C0xz + C2xz)/30.;
  Gek_y[7][0] = -(C0xz + 2.*C3xz)/30.;
  Gek_y[7][1] = -(C0xz + C3xz)/30.;
  Gek_y[7][2] = -(C0xz + C3xz)/30.;
  Gek_y[7][3] = -(2.*C0xz + C3xz)/30.;
  Gek_y[8][0] = -(C1xz + C3xz)/30.;
  Gek_y[8][1] = -(C1xz + 2.*C3xz)/30.;
  Gek_y[8][2] = -(C1xz + C3xz)/30.;
  Gek_y[8][3] = -(2.*C1xz + C3xz)/30.;
  Gek_y[9][0] = -(C2xz + C3xz)/30.;
  Gek_y[9][1] = -(C2xz + C3xz)/30.;
  Gek_y[9][2] = -(C2xz + 2.*C3xz)/30.;
  Gek_y[9][3] = -(2.*C2xz + C3xz)/30.;

  Gek_z[0][0] = -C0xy/40.;
  Gek_z[0][1] = C0xy/120.;
  Gek_z[0][2] = C0xy/120.;
  Gek_z[0][3] = C0xy/120.;
  Gek_z[1][0] = C1xy/120.;
  Gek_z[1][1] = -C1xy/40.;
  Gek_z[1][2] = C1xy/120.;
  Gek_z[1][3] = C1xy/120.;
  Gek_z[2][0] = C2xy/120.;
  Gek_z[2][1] = C2xy/120.;
  Gek_z[2][2] = -C2xy/40.;
  Gek_z[2][3] = C2xy/120.;
  Gek_z[3][0] = C3xy/120.;
  Gek_z[3][1] = C3xy/120.;
  Gek_z[3][2] = C3xy/120.;
  Gek_z[3][3] = -C3xy/40.;
  Gek_z[4][0] = -(C0xy + 2.*C1xy)/30.;
  Gek_z[4][1] = -(2.*C0xy + C1xy)/30.;
  Gek_z[4][2] = -(C0xy + C1xy)/30.;
  Gek_z[4][3] = -(C0xy + C1xy)/30.;
  Gek_z[5][0] = -(C1xy + C2xy)/30.;
  Gek_z[5][1] = -(C1xy + 2.*C2xy)/30.;
  Gek_z[5][2] = -(2.*C1xy + C2xy)/30.;
  Gek_z[5][3] = -(C1xy + C2xy)/30.;
  Gek_z[6][0] = -(C0xy + 2.*C2xy)/30.;
  Gek_z[6][1] = -(C0xy + C2xy)/30.;
  Gek_z[6][2] = -(2.*C0xy + C2xy)/30.;
  Gek_z[6][3] = -(C0xy + C2xy)/30.;
  Gek_z[7][0] = -(C0xy + 2.*C3xy)/30.;
  Gek_z[7][1] = -(C0xy + C3xy)/30.;
  Gek_z[7][2] = -(C0xy + C3xy)/30.;
  Gek_z[7][3] = -(2.*C0xy + C3xy)/30.;
  Gek_z[8][0] = -(C1xy + C3xy)/30.;
  Gek_z[8][1] = -(C1xy + 2.*C3xy)/30.;
  Gek_z[8][2] = -(C1xy + C3xy)/30.;
  Gek_z[8][3] = -(2.*C1xy + C3xy)/30.;
  Gek_z[9][0] = -(C2xy + C3xy)/30.;
  Gek_z[9][1] = -(C2xy + C3xy)/30.;
  Gek_z[9][2] = -(C2xy + 2.*C3xy)/30.;
  Gek_z[9][3] = -(2.*C2xy + C3xy)/30.;

	return 0;
}



PetscErrorCode Calc3DP1P2_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscInt n3 = e_inter[en][3];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];

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

  Gek_x[0][0] = -C0yz/40.;
  Gek_x[1][0] =  C0yz/120.;
  Gek_x[2][0] =  C0yz/120.;
  Gek_x[3][0] =  C0yz/120.;
  Gek_x[0][1] =  C1yz/120.;
  Gek_x[1][1] = -C1yz/40.;
  Gek_x[2][1] =  C1yz/120.;
  Gek_x[3][1] =  C1yz/120.;
  Gek_x[0][2] =  C2yz/120.;
  Gek_x[1][2] =  C2yz/120.;
  Gek_x[2][2] = -C2yz/40.;
  Gek_x[3][2] =  C2yz/120.;
  Gek_x[0][3] =  C3yz/120.;
  Gek_x[1][3] =  C3yz/120.;
  Gek_x[2][3] =  C3yz/120.;
  Gek_x[3][3] = -C3yz/40.;
  Gek_x[0][4] = -(C0yz + 2.*C1yz)/30.;
  Gek_x[1][4] = -(2.*C0yz + C1yz)/30.;
  Gek_x[2][4] = -(C0yz + C1yz)/30.;
  Gek_x[3][4] = -(C0yz + C1yz)/30.;
  Gek_x[0][5] = -(C1yz + C2yz)/30.;
  Gek_x[1][5] = -(C1yz + 2.*C2yz)/30.;
  Gek_x[2][5] = -(2.*C1yz + C2yz)/30.;
  Gek_x[3][5] = -(C1yz + C2yz)/30.;
  Gek_x[0][6] = -(C0yz + 2.*C2yz)/30.;
  Gek_x[1][6] = -(C0yz + C2yz)/30.;
  Gek_x[2][6] = -(2.*C0yz + C2yz)/30.;
  Gek_x[3][6] = -(C0yz + C2yz)/30.;
  Gek_x[0][7] = -(C0yz + 2.*C3yz)/30.;
  Gek_x[1][7] = -(C0yz + C3yz)/30.;
  Gek_x[2][7] = -(C0yz + C3yz)/30.;
  Gek_x[3][7] = -(2.*C0yz + C3yz)/30.;
  Gek_x[0][8] = -(C1yz + C3yz)/30.;
  Gek_x[1][8] = -(C1yz + 2.*C3yz)/30.;
  Gek_x[2][8] = -(C1yz + C3yz)/30.;
  Gek_x[3][8] = -(2.*C1yz + C3yz)/30.;
  Gek_x[0][9] = -(C2yz + C3yz)/30.;
  Gek_x[1][9] = -(C2yz + C3yz)/30.;
  Gek_x[2][9] = -(C2yz + 2.*C3yz)/30.;
  Gek_x[3][9] = -(2.*C2yz + C3yz)/30.;

  Gek_y[0][0] = -C0xz/40.;
  Gek_y[1][0] = C0xz/120.;
  Gek_y[2][0] = C0xz/120.;
  Gek_y[3][0] = C0xz/120.;
  Gek_y[0][1] = C1xz/120.;
  Gek_y[1][1] = -C1xz/40.;
  Gek_y[2][1] = C1xz/120.;
  Gek_y[3][1] = C1xz/120.;
  Gek_y[0][2] = C2xz/120.;
  Gek_y[1][2] = C2xz/120.;
  Gek_y[2][2] = -C2xz/40.;
  Gek_y[3][2] = C2xz/120.;
  Gek_y[0][3] = C3xz/120.;
  Gek_y[1][3] = C3xz/120.;
  Gek_y[2][3] = C3xz/120.;
  Gek_y[3][3] = -C3xz/40.;
  Gek_y[0][4] = -(C0xz + 2.*C1xz)/30.;
  Gek_y[1][4] = -(2.*C0xz + C1xz)/30.;
  Gek_y[2][4] = -(C0xz + C1xz)/30.;
  Gek_y[3][4] = -(C0xz + C1xz)/30.;
  Gek_y[0][5] = -(C1xz + C2xz)/30.;
  Gek_y[1][5] = -(C1xz + 2.*C2xz)/30.;
  Gek_y[2][5] = -(2.*C1xz + C2xz)/30.;
  Gek_y[3][5] = -(C1xz + C2xz)/30.;
  Gek_y[0][6] = -(C0xz + 2.*C2xz)/30.;
  Gek_y[1][6] = -(C0xz + C2xz)/30.;
  Gek_y[2][6] = -(2.*C0xz + C2xz)/30.;
  Gek_y[3][6] = -(C0xz + C2xz)/30.;
  Gek_y[0][7] = -(C0xz + 2.*C3xz)/30.;
  Gek_y[1][7] = -(C0xz + C3xz)/30.;
  Gek_y[2][7] = -(C0xz + C3xz)/30.;
  Gek_y[3][7] = -(2.*C0xz + C3xz)/30.;
  Gek_y[0][8] = -(C1xz + C3xz)/30.;
  Gek_y[1][8] = -(C1xz + 2.*C3xz)/30.;
  Gek_y[2][8] = -(C1xz + C3xz)/30.;
  Gek_y[3][8] = -(2.*C1xz + C3xz)/30.;
  Gek_y[0][9] = -(C2xz + C3xz)/30.;
  Gek_y[1][9] = -(C2xz + C3xz)/30.;
  Gek_y[2][9] = -(C2xz + 2.*C3xz)/30.;
  Gek_y[3][9] = -(2.*C2xz + C3xz)/30.;

  Gek_z[0][0] = -C0xy/40.;
  Gek_z[1][0] = C0xy/120.;
  Gek_z[2][0] = C0xy/120.;
  Gek_z[3][0] = C0xy/120.;
  Gek_z[0][1] = C1xy/120.;
  Gek_z[1][1] = -C1xy/40.;
  Gek_z[2][1] = C1xy/120.;
  Gek_z[3][1] = C1xy/120.;
  Gek_z[0][2] = C2xy/120.;
  Gek_z[1][2] = C2xy/120.;
  Gek_z[2][2] = -C2xy/40.;
  Gek_z[3][2] = C2xy/120.;
  Gek_z[0][3] = C3xy/120.;
  Gek_z[1][3] = C3xy/120.;
  Gek_z[2][3] = C3xy/120.;
  Gek_z[3][3] = -C3xy/40.;
  Gek_z[0][4] = -(C0xy + 2.*C1xy)/30.;
  Gek_z[1][4] = -(2.*C0xy + C1xy)/30.;
  Gek_z[2][4] = -(C0xy + C1xy)/30.;
  Gek_z[3][4] = -(C0xy + C1xy)/30.;
  Gek_z[0][5] = -(C1xy + C2xy)/30.;
  Gek_z[1][5] = -(C1xy + 2.*C2xy)/30.;
  Gek_z[2][5] = -(2.*C1xy + C2xy)/30.;
  Gek_z[3][5] = -(C1xy + C2xy)/30.;
  Gek_z[0][6] = -(C0xy + 2.*C2xy)/30.;
  Gek_z[1][6] = -(C0xy + C2xy)/30.;
  Gek_z[2][6] = -(2.*C0xy + C2xy)/30.;
  Gek_z[3][6] = -(C0xy + C2xy)/30.;
  Gek_z[0][7] = -(C0xy + 2.*C3xy)/30.;
  Gek_z[1][7] = -(C0xy + C3xy)/30.;
  Gek_z[2][7] = -(C0xy + C3xy)/30.;
  Gek_z[3][7] = -(2.*C0xy + C3xy)/30.;
  Gek_z[0][8] = -(C1xy + C3xy)/30.;
  Gek_z[1][8] = -(C1xy + 2.*C3xy)/30.;
  Gek_z[2][8] = -(C1xy + C3xy)/30.;
  Gek_z[3][8] = -(2.*C1xy + C3xy)/30.;
  Gek_z[0][9] = -(C2xy + C3xy)/30.;
  Gek_z[1][9] = -(C2xy + C3xy)/30.;
  Gek_z[2][9] = -(C2xy + 2.*C3xy)/30.;
  Gek_z[3][9] = -(2.*C2xy + C3xy)/30.;

	return 0;
}

PetscErrorCode Calc3DP1bP1_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscInt n3 = e_inter[en][3];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];

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

  Gek_x[0][0] =  -137.*C0yz/2520; 
  Gek_x[0][1] =  -(73.*C0yz - 32.*C2yz - 32.*C3yz)/2520;
  Gek_x[0][2] =  -(73.*C0yz - 32.*C1yz - 32.*C3yz)/2520;
  Gek_x[0][3] =  -(73.*C0yz - 32.*C1yz - 32.*C2yz)/2520;

  Gek_x[1][0] =  -(73.*C1yz - 32.*C2yz - 32.*C3yz)/2520;
  Gek_x[1][1] =  -137.*C1yz/2520;
  Gek_x[1][2] =  (32.*C0yz - 73.*C1yz + 32.*C3yz)/2520;
  Gek_x[1][3] =  (32.*C0yz - 73.*C1yz + 32.*C2yz)/2520;

  Gek_x[2][0] =  (32.*C1yz - 73.*C2yz + 32.*C3yz)/2520;
  Gek_x[2][1] =  (32.*C0yz - 73.*C2yz + 32.*C3yz)/2520;
  Gek_x[2][2] =  -137.*C2yz/2520;
  Gek_x[2][3] =  (32.*C0yz + 32.*C1yz - 73.*C2yz)/2520;

  Gek_x[3][0] =  (32.*C1yz + 32.*C2yz - 73.*C3yz)/2520;
  Gek_x[3][1] =  (32.*C0yz + 32.*C2yz - 73.*C3yz)/2520;
  Gek_x[3][2] =  (32.*C0yz + 32.*C1yz - 73.*C3yz)/2520;
  Gek_x[3][3] =  -137.*C3yz/2520;

  Gek_x[4][0] =  16.*C0yz/315;
  Gek_x[4][1] =  16.*C1yz/315;
  Gek_x[4][2] =  16.*C2yz/315;
  Gek_x[4][3] =  16.*C3yz/315;
  

  Gek_y[0][0] =  -137.*C0xz/2520;
  Gek_y[0][1] =  -(73.*C0xz - 32.*C2xz - 32.*C3xz)/2520;
  Gek_y[0][2] =  -(73.*C0xz - 32.*C1xz - 32.*C3xz)/2520;
  Gek_y[0][3] =  -(73.*C0xz - 32.*C1xz - 32.*C2xz)/2520;
  
  Gek_y[1][0] =  -(73.*C1xz - 32.*C2xz - 32.*C3xz)/2520;
  Gek_y[1][1] =  -137.*C1xz/2520;
  Gek_y[1][2] =  (32.*C0xz - 73.*C1xz + 32.*C3xz)/2520;
  Gek_y[1][3] =  (32.*C0xz - 73.*C1xz + 32.*C2xz)/2520;
  
  Gek_y[2][0] =  (32.*C1xz - 73.*C2xz + 32.*C3xz)/2520;
  Gek_y[2][1] =  (32.*C0xz - 73.*C2xz + 32.*C3xz)/2520;
  Gek_y[2][2] =  -137.*C2xz/2520;
  Gek_y[2][3] =  (32.*C0xz + 32.*C1xz - 73.*C2xz)/2520;

  Gek_y[3][0] =  (32.*C1xz + 32.*C2xz - 73.*C3xz)/2520;
  Gek_y[3][1] =  (32.*C0xz + 32.*C2xz - 73.*C3xz)/2520;
  Gek_y[3][2] =  (32.*C0xz + 32.*C1xz - 73.*C3xz)/2520;
  Gek_y[3][3] =  -137.*C3xz/2520;
  
  Gek_y[4][0] =  16.*C0xz/315;
  Gek_y[4][1] =  16.*C1xz/315;
  Gek_y[4][2] =  16.*C2xz/315;
  Gek_y[4][3] =  16.*C3xz/315;
  
  
  Gek_z[0][0] =  -137.*C0xy/2520;
  Gek_z[0][1] =  -(73.*C0xy - 32.*C2xy - 32.*C3xy)/2520;
  Gek_z[0][2] =  -(73.*C0xy - 32.*C1xy - 32.*C3xy)/2520;
  Gek_z[0][3] =  -(73.*C0xy - 32.*C1xy - 32.*C2xy)/2520;

  Gek_z[1][0] =  -(73.*C1xy - 32.*C2xy - 32.*C3xy)/2520;
  Gek_z[1][1] =  -137.*C1xy/2520;
  Gek_z[1][2] =  (32.*C0xy - 73.*C1xy + 32.*C3xy)/2520;
  Gek_z[1][3] =  (32.*C0xy - 73.*C1xy + 32.*C2xy)/2520;

  Gek_z[2][0] =  (32.*C1xy - 73.*C2xy + 32.*C3xy)/2520;
  Gek_z[2][1] =  (32.*C0xy - 73.*C2xy + 32.*C3xy)/2520;
  Gek_z[2][2] =  -137.*C2xy/2520;
  Gek_z[2][3] =  (32.*C0xy + 32.*C1xy - 73.*C2xy)/2520;

  Gek_z[3][0] =  (32.*C1xy + 32.*C2xy - 73.*C3xy)/2520;
  Gek_z[3][1] =  (32.*C0xy + 32.*C2xy - 73.*C3xy)/2520;
  Gek_z[3][2] =  (32.*C0xy + 32.*C1xy - 73.*C3xy)/2520;
  Gek_z[3][3] =  -137.*C3xy/2520;
  
  Gek_z[4][0] =  16.*C0xy/315;
  Gek_z[4][1] =  16.*C1xy/315;
  Gek_z[4][2] =  16.*C2xy/315;
  Gek_z[4][3] =  16.*C3xy/315;

  return 0;
}


PetscErrorCode Calc3DP1P1b_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscInt n3 = e_inter[en][3];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];

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

  Gek_x[0][0] =  -137.*C0yz/2520; 
  Gek_x[1][0] =  -(73.*C0yz - 32.*C2yz - 32.*C3yz)/2520;
  Gek_x[2][0] =  -(73.*C0yz - 32.*C1yz - 32.*C3yz)/2520;
  Gek_x[3][0] =  -(73.*C0yz - 32.*C1yz - 32.*C2yz)/2520;
  Gek_x[0][1] =  -(73.*C1yz - 32.*C2yz - 32.*C3yz)/2520;
  Gek_x[1][1] =  -137.*C1yz/2520;
  Gek_x[2][1] =  (32.*C0yz - 73.*C1yz + 32.*C3yz)/2520;
  Gek_x[3][1] =  (32.*C0yz - 73.*C1yz + 32.*C2yz)/2520;
  Gek_x[0][2] =  (32.*C1yz - 73.*C2yz + 32.*C3yz)/2520;
  Gek_x[1][2] =  (32.*C0yz - 73.*C2yz + 32.*C3yz)/2520;
  Gek_x[2][2] =  -137.*C2yz/2520;
  Gek_x[3][2] =  (32.*C0yz + 32.*C1yz - 73.*C2yz)/2520;
  Gek_x[0][3] =  (32.*C1yz + 32.*C2yz - 73.*C3yz)/2520;
  Gek_x[1][3] =  (32.*C0yz + 32.*C2yz - 73.*C3yz)/2520;
  Gek_x[2][3] =  (32.*C0yz + 32.*C1yz - 73.*C3yz)/2520;
  Gek_x[3][3] =  -137.*C3yz/2520;
  Gek_x[0][4] =  16.*C0yz/315;
  Gek_x[1][4] =  16.*C1yz/315;
  Gek_x[2][4] =  16.*C2yz/315;
  Gek_x[3][4] =  16.*C3yz/315;
  
  Gek_y[0][0] =  -137.*C0xz/2520;
  Gek_y[1][0] =  -(73.*C0xz - 32.*C2xz - 32.*C3xz)/2520;
  Gek_y[2][0] =  -(73.*C0xz - 32.*C1xz - 32.*C3xz)/2520;
  Gek_y[3][0] =  -(73.*C0xz - 32.*C1xz - 32.*C2xz)/2520;
  Gek_y[0][1] =  -(73.*C1xz - 32.*C2xz - 32.*C3xz)/2520;
  Gek_y[1][1] =  -137.*C1xz/2520;
  Gek_y[2][1] =  (32.*C0xz - 73.*C1xz + 32.*C3xz)/2520;
  Gek_y[3][1] =  (32.*C0xz - 73.*C1xz + 32.*C2xz)/2520;
  Gek_y[0][2] =  (32.*C1xz - 73.*C2xz + 32.*C3xz)/2520;
  Gek_y[1][2] =  (32.*C0xz - 73.*C2xz + 32.*C3xz)/2520;
  Gek_y[2][2] =  -137.*C2xz/2520;
  Gek_y[3][2] =  (32.*C0xz + 32.*C1xz - 73.*C2xz)/2520;
  Gek_y[0][3] =  (32.*C1xz + 32.*C2xz - 73.*C3xz)/2520;
  Gek_y[1][3] =  (32.*C0xz + 32.*C2xz - 73.*C3xz)/2520;
  Gek_y[2][3] =  (32.*C0xz + 32.*C1xz - 73.*C3xz)/2520;
  Gek_y[3][3] =  -137.*C3xz/2520;
  Gek_y[0][4] =  16.*C0xz/315;
  Gek_y[1][4] =  16.*C1xz/315;
  Gek_y[2][4] =  16.*C2xz/315;
  Gek_y[3][4] =  16.*C3xz/315;
  
  Gek_z[0][0] =  -137.*C0xy/2520;
  Gek_z[1][0] =  -(73.*C0xy - 32.*C2xy - 32.*C3xy)/2520;
  Gek_z[2][0] =  -(73.*C0xy - 32.*C1xy - 32.*C3xy)/2520;
  Gek_z[3][0] =  -(73.*C0xy - 32.*C1xy - 32.*C2xy)/2520;
  Gek_z[0][1] =  -(73.*C1xy - 32.*C2xy - 32.*C3xy)/2520;
  Gek_z[1][1] =  -137.*C1xy/2520;
  Gek_z[2][1] =  (32.*C0xy - 73.*C1xy + 32.*C3xy)/2520;
  Gek_z[3][1] =  (32.*C0xy - 73.*C1xy + 32.*C2xy)/2520;
  Gek_z[0][2] =  (32.*C1xy - 73.*C2xy + 32.*C3xy)/2520;
  Gek_z[1][2] =  (32.*C0xy - 73.*C2xy + 32.*C3xy)/2520;
  Gek_z[2][2] =  -137.*C2xy/2520;
  Gek_z[3][2] =  (32.*C0xy + 32.*C1xy - 73.*C2xy)/2520;
  Gek_z[0][3] =  (32.*C1xy + 32.*C2xy - 73.*C3xy)/2520;
  Gek_z[1][3] =  (32.*C0xy + 32.*C2xy - 73.*C3xy)/2520;
  Gek_z[2][3] =  (32.*C0xy + 32.*C1xy - 73.*C3xy)/2520;
  Gek_z[3][3] =  -137.*C3xy/2520;
  Gek_z[0][4] =  16.*C0xy/315;
  Gek_z[1][4] =  16.*C1xy/315;
  Gek_z[2][4] =  16.*C2xy/315;
  Gek_z[3][4] =  16.*C3xy/315;

  return 0;
}



PetscErrorCode Calc3DP1bP1b_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  
  const PetscInt n0 = e_inter[en][0];
  const PetscInt n1 = e_inter[en][1];
  const PetscInt n2 = e_inter[en][2];
  const PetscInt n3 = e_inter[en][3];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];

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

  Gek_x[0][0] =  -C0yz/24.;
  Gek_x[0][1] =  -(73.*C0yz + 32.*C1yz)/2520.;
  Gek_x[0][2] =  -(73.*C0yz + 32.*C2yz)/2520.;
  Gek_x[0][3] =  -(73.*C0yz + 32.*C3yz)/2520.;
  Gek_x[0][4] =  -16.*C0yz/315.;
  Gek_x[1][0] =  -(32.*C0yz + 73.*C1yz)/2520.;
  Gek_x[1][1] =  -C1yz/24.;
  Gek_x[1][2] =  -(73.*C1yz + 32.*C2yz)/2520.;
  Gek_x[1][3] =  -(73.*C1yz + 32.*C3yz)/2520.;
  Gek_x[1][4] =  -16.*C1yz/315.;
  Gek_x[2][0] =  -(32.*C0yz + 73.*C2yz)/2520.;
  Gek_x[2][1] =  -(32.*C1yz + 73.*C2yz)/2520.;
  Gek_x[2][2] =  -C2yz/24.;
  Gek_x[2][3] =  -(73.*C2yz + 32.*C3yz)/2520.;
  Gek_x[2][4] =  -16.*C2yz/315.;
  Gek_x[3][0] =  -(32.*C0yz + 73.*C3yz)/2520.;
  Gek_x[3][1] =  -(32.*C1yz + 73.*C3yz)/2520.;
  Gek_x[3][2] =  -(32.*C2yz + 73.*C3yz)/2520.;
  Gek_x[3][3] =  -C3yz/24.;
  Gek_x[3][4] =  -16.*C3yz/315.;
  Gek_x[4][0] =   16.*C0yz/315.;
  Gek_x[4][1] =   16.*C1yz/315.;
  Gek_x[4][2] =   16.*C2yz/315.;
  Gek_x[4][3] =   16.*C3yz/315.;
  Gek_x[4][4] =   0.;

  Gek_y[0][0] =  -C0xz/24.;
  Gek_y[0][1] =  -(73.*C0xz + 32.*C1xz)/2520.;
  Gek_y[0][2] =  -(73.*C0xz + 32.*C2xz)/2520.;
  Gek_y[0][3] =  -(73.*C0xz + 32.*C3xz)/2520.;
  Gek_y[0][4] =  -16.*C0xz/315.;
  Gek_y[1][0] =  -(32.*C0xz + 73.*C1xz)/2520.;
  Gek_y[1][1] =  -C1xz/24.;
  Gek_y[1][2] =  -(73.*C1xz + 32.*C2xz)/2520.;
  Gek_y[1][3] =  -(73.*C1xz + 32.*C3xz)/2520.;
  Gek_y[1][4] =  -16.*C1xz/315.;
  Gek_y[2][0] =  -(32.*C0xz + 73.*C2xz)/2520.;
  Gek_y[2][1] =  -(32.*C1xz + 73.*C2xz)/2520.;
  Gek_y[2][2] =  -C2xz/24.;
  Gek_y[2][3] =  -(73.*C2xz + 32.*C3xz)/2520.;
  Gek_y[2][4] =  -16.*C2xz/315.;
  Gek_y[3][0] =  -(32.*C0xz + 73.*C3xz)/2520.;
  Gek_y[3][1] =  -(32.*C1xz + 73.*C3xz)/2520.;
  Gek_y[3][2] =  -(32.*C2xz + 73.*C3xz)/2520.;
  Gek_y[3][3] =  -C3xz/24.;
  Gek_y[3][4] =  -16.*C3xz/315.;
  Gek_y[4][0] =   16.*C0xz/315.;
  Gek_y[4][1] =   16.*C1xz/315.;
  Gek_y[4][2] =   16.*C2xz/315.;
  Gek_y[4][3] =   16.*C3xz/315.;
  Gek_y[4][4] =   0.;

  Gek_z[0][0] =  -C0xy/24.;
  Gek_z[0][1] =  -(73.*C0xy + 32.*C1xy)/2520.;
  Gek_z[0][2] =  -(73.*C0xy + 32.*C2xy)/2520.;
  Gek_z[0][3] =  -(73.*C0xy + 32.*C3xy)/2520.;
  Gek_z[0][4] =  -16.*C0xy/315.;
  Gek_z[1][0] =  -(32.*C0xy + 73.*C1xy)/2520.;
  Gek_z[1][1] =  -C1xy/24.;
  Gek_z[1][2] =  -(73.*C1xy + 32.*C2xy)/2520.;
  Gek_z[1][3] =  -(73.*C1xy + 32.*C3xy)/2520.;
  Gek_z[1][4] =  -16.*C1xy/315.;
  Gek_z[2][0] =  -(32.*C0xy + 73.*C2xy)/2520.;
  Gek_z[2][1] =  -(32.*C1xy + 73.*C2xy)/2520.;
  Gek_z[2][2] =  -C2xy/24.;
  Gek_z[2][3] =  -(73.*C2xy + 32.*C3xy)/2520.;
  Gek_z[2][4] =  -16.*C2xy/315.;
  Gek_z[3][0] =  -(32.*C0xy + 73.*C3xy)/2520.;
  Gek_z[3][1] =  -(32.*C1xy + 73.*C3xy)/2520.;
  Gek_z[3][2] =  -(32.*C2xy + 73.*C3xy)/2520.;
  Gek_z[3][3] =  -C3xy/24.;
  Gek_z[3][4] =  -16.*C3xy/315.;
  Gek_z[4][0] =   16.*C0xy/315.;
  Gek_z[4][1] =   16.*C1xy/315.;
  Gek_z[4][2] =   16.*C2xy/315.;
  Gek_z[4][3] =   16.*C3xy/315.;
  Gek_z[4][4] =   0.;

  return 0;
}

#endif