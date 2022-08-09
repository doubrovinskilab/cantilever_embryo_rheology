#ifndef ElementalMatrixDiffusion_C
#define ElementalMatrixDiffusion_C

PetscErrorCode Calc2DP1_Dek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Dek)
{

	const PetscInt n0 = e_inter[en][0];
	const PetscInt n1 = e_inter[en][1];
	const PetscInt n2 = e_inter[en][2];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2];
	
	const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;
	const PetscScalar Inv_2Del = (1.0/(2.0*Delta));
	
	Dek[0][0] = Inv_2Del* (  x1*x1 - 2.*x1*x2 + x2*x2 + y1*y1 - 2.*y1*y2 + y2*y2);
	Dek[0][1] = Inv_2Del* (-(x0*x1 - x0*x2 - x1*x2 + x2*x2 + y0*y1 - y0*y2 - y1*y2 + y2*y2));
	Dek[0][2] = Inv_2Del* (  x0*x1 - x0*x2 - x1*x1 + x1*x2 + y0*y1 - y0*y2 - y1*y1 + y1*y2);
	Dek[1][0] = Inv_2Del* (-(x0*x1 - x0*x2 - x1*x2 + x2*x2 + y0*y1 - y0*y2 - y1*y2 + y2*y2));
	Dek[1][1] = Inv_2Del* (  x0*x0 - 2.*x0*x2 + x2*x2 + y0*y0 - 2.*y0*y2 + y2*y2);
	Dek[1][2] = Inv_2Del* (-(x0*x0 - x0*x1 - x0*x2 + x1*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y2));
	Dek[2][0] = Inv_2Del* (  x0*x1 - x0*x2 - x1*x1 + x1*x2 + y0*y1 - y0*y2 - y1*y1 + y1*y2);
	Dek[2][1] = Inv_2Del* (-(x0*x0 - x0*x1 - x0*x2 + x1*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y2));
	Dek[2][2] = Inv_2Del* (  x0*x0 - 2.*x0*x1 + x1*x1 + y0*y0 - 2.*y0*y1 + y1*y1);

	return 0;
}

PetscErrorCode Calc2DP2_Dek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **eP2_inter, PetscScalar **Dek)
{

	const PetscInt n0 = eP2_inter[en][0];
	const PetscInt n1 = eP2_inter[en][1];
	const PetscInt n2 = eP2_inter[en][2];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2];
	
	const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;

	Dek[0][0] = (x1*x1 - 2.*x1*x2 + x2*x2 + y1*y1 - 2.*y1*y2 + y2*y2)/(2.*Delta);
	Dek[0][1] = (x0*x1 - x0*x2 - x1*x2 + x2*x2 + y0*y1 - y0*y2 - y1*y2 + y2*y2)/(6.*Delta);
	Dek[0][2] = -(x0*x1 - x0*x2 - x1*x1 + x1*x2 + y0*y1 - y0*y2 - y1*y1 + y1*y2)/(6.*Delta);
	Dek[0][3] = -2.*(x0*x1 - x0*x2 - x1*x2 + x2*x2 + y0*y1 - y0*y2 - y1*y2 + y2*y2)/(3.*Delta);
	Dek[0][4] = 0.;
	Dek[0][5] = 2.*(x0*x1 - x0*x2 - x1*x1 + x1*x2 + y0*y1 - y0*y2 - y1*y1 + y1*y2)/(3.*Delta);
	Dek[1][0] = (x0*x1 - x0*x2 - x1*x2 + x2*x2 + y0*y1 - y0*y2 - y1*y2 + y2*y2)/(6.*Delta);
	Dek[1][1] = (x0*x0 - 2.*x0*x2 + x2*x2 + y0*y0 - 2.*y0*y2 + y2*y2)/(2.*Delta);
	Dek[1][2] = (x0*x0 - x0*x1 - x0*x2 + x1*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y2)/(6.*Delta);
	Dek[1][3] = -2.*(x0*x1 - x0*x2 - x1*x2 + x2*x2 + y0*y1 - y0*y2 - y1*y2 + y2*y2)/(3.*Delta);
	Dek[1][4] = -2.*(x0*x0 - x0*x1 - x0*x2 + x1*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y2)/(3.*Delta);
	Dek[1][5] = 0.;
	Dek[2][0] = -(x0*x1 - x0*x2 - x1*x1 + x1*x2 + y0*y1 - y0*y2 - y1*y1 + y1*y2)/(6.*Delta);
	Dek[2][1] = (x0*x0 - x0*x1 - x0*x2 + x1*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y2)/(6.*Delta);
	Dek[2][2] = (x0*x0 - 2.*x0*x1 + x1*x1 + y0*y0 - 2.*y0*y1 + y1*y1)/(2.*Delta);
	Dek[2][3] = 0.;
	Dek[2][4] = -2.*(x0*x0 - x0*x1 - x0*x2 + x1*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y2)/(3.*Delta);
	Dek[2][5] = 2.*(x0*x1 - x0*x2 - x1*x1 + x1*x2 + y0*y1 - y0*y2 - y1*y1 + y1*y2)/(3.*Delta);
	Dek[3][0] = -2.*(x0*x1 - x0*x2 - x1*x2 + x2*x2 + y0*y1 - y0*y2 - y1*y2 + y2*y2)/(3.*Delta);
	Dek[3][1] = -2.*(x0*x1 - x0*x2 - x1*x2 + x2*x2 + y0*y1 - y0*y2 - y1*y2 + y2*y2)/(3.*Delta);
	Dek[3][2] = 0.;
	Dek[3][3] = 4.*(x0*x0 - x0*x1 - x0*x2 + x1*x1 - x1*x2 + x2*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y1 - y1*y2 + y2*y2)/(3.*Delta);
	Dek[3][4] = 4.*(x0*x1 - x0*x2 - x1*x1 + x1*x2 + y0*y1 - y0*y2 - y1*y1 + y1*y2)/(3.*Delta);
	Dek[3][5] = -4.*(x0*x0 - x0*x1 - x0*x2 + x1*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y2)/(3.*Delta);
	Dek[4][0] = 0.;
	Dek[4][1] = -2.*(x0*x0 - x0*x1 - x0*x2 + x1*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y2)/(3.*Delta);
	Dek[4][2] = -2.*(x0*x0 - x0*x1 - x0*x2 + x1*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y2)/(3.*Delta);
	Dek[4][3] = 4.*(x0*x1 - x0*x2 - x1*x1 + x1*x2 + y0*y1 - y0*y2 - y1*y1 + y1*y2)/(3.*Delta);
	Dek[4][4] = 4.*(x0*x0 - x0*x1 - x0*x2 + x1*x1 - x1*x2 + x2*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y1 - y1*y2 + y2*y2)/(3.*Delta);
	Dek[4][5] = -4.*(x0*x1 - x0*x2 - x1*x2 + x2*x2 + y0*y1 - y0*y2 - y1*y2 + y2*y2)/(3.*Delta);
	Dek[5][0] = 2.*(x0*x1 - x0*x2 - x1*x1 + x1*x2 + y0*y1 - y0*y2 - y1*y1 + y1*y2)/(3.*Delta);
	Dek[5][1] = 0.;
	Dek[5][2] = 2.*(x0*x1 - x0*x2 - x1*x1 + x1*x2 + y0*y1 - y0*y2 - y1*y1 + y1*y2)/(3.*Delta);
	Dek[5][3] = -4.*(x0*x0 - x0*x1 - x0*x2 + x1*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y2)/(3.*Delta);
	Dek[5][4] = -4.*(x0*x1 - x0*x2 - x1*x2 + x2*x2 + y0*y1 - y0*y2 - y1*y2 + y2*y2)/(3.*Delta);
	Dek[5][5] =  4.*(x0*x0 - x0*x1 - x0*x2 + x1*x1 - x1*x2 + x2*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y1 - y1*y2 + y2*y2)/(3.*Delta);

	return 0;
}

PetscErrorCode Calc2DP1b_Dek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **eP2_inter, PetscScalar **Dek)
{

	const PetscInt n0 = eP2_inter[en][0];
	const PetscInt n1 = eP2_inter[en][1];
	const PetscInt n2 = eP2_inter[en][2];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2];
	
	const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;

  Dek[0][0] = (9.*x0*x0 - 9.*x0*x1 - 9.*x0*x2 + 19.*x1*x1 - 29.*x1*x2 + 19.*x2*x2 + 9.*y0*y0 - 9.*y0*y1 - 9.*y0*y2 + 19.*y1*y1 - 29.*y1*y2 + 19.*y2*y2)/(20.*Delta);
  Dek[0][1] = (9.*x0*x0 - 19.*x0*x1 + x0*x2 + 9.*x1*x1 + x1*x2 - x2*x2 + 9.*y0*y0 - 19.*y0*y1 + y0*y2 + 9.*y1*y1 + y1*y2 - y2*y2)/(20.*Delta);
  Dek[0][2] = (9.*x0*x0 + x0*x1 - 19.*x0*x2 - x1*x1 + x1*x2 + 9.*x2*x2 + 9.*y0*y0 + y0*y1 - 19.*y0*y2 - y1*y1 + y1*y2 + 9.*y2*y2)/(20.*Delta);
  Dek[0][3] = -27.*(x0*x0 - x0*x1 - x0*x2 + x1*x1 - x1*x2 + x2*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y1 - y1*y2 + y2*y2)/(20.*Delta);

  Dek[1][0] = (9.*x0*x0 - 19.*x0*x1 + x0*x2 + 9.*x1*x1 + x1*x2 - x2*x2 + 9.*y0*y0 - 19.*y0*y1 + y0*y2 + 9.*y1*y1 + y1*y2 - y2*y2)/(20.*Delta);
  Dek[1][1] =  (19.*x0*x0 - 9.*x0*x1 - 29.*x0*x2 + 9.*x1*x1 - 9.*x1*x2 + 19.*x2*x2 + 19.*y0*y0 - 9.*y0*y1 - 29.*y0*y2 + 9.*y1*y1 - 9.*y1*y2 + 19.*y2*y2)/(20.*Delta);
  Dek[1][2] = -(x0*x0 - x0*x1 - x0*x2 - 9.*x1*x1 + 19.*x1*x2 - 9.*x2*x2 + y0*y0 - y0*y1 - y0*y2 - 9.*y1*y1 + 19.*y1*y2 - 9.*y2*y2)/(20.*Delta);
  Dek[1][3] = -27.*(x0*x0 - x0*x1 - x0*x2 + x1*x1 - x1*x2 + x2*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y1 - y1*y2 + y2*y2)/(20.*Delta);

  Dek[2][0] = (9.*x0*x0 + x0*x1 - 19.*x0*x2 - x1*x1 + x1*x2 + 9.*x2*x2 + 9.*y0*y0 + y0*y1 - 19.*y0*y2 - y1*y1 + y1*y2 + 9.*y2*y2)/(20.*Delta);
  Dek[2][1] = -(x0*x0 - x0*x1 - x0*x2 - 9.*x1*x1 + 19.*x1*x2 - 9.*x2*x2 + y0*y0 - y0*y1 - y0*y2 - 9.*y1*y1 + 19.*y1*y2 - 9.*y2*y2)/(20.*Delta);
  Dek[2][2] = (19.*x0*x0 - 29.*x0*x1 - 9.*x0*x2 + 19.*x1*x1 - 9.*x1*x2 + 9.*x2*x2 + 19.*y0*y0 - 29.*y0*y1 - 9.*y0*y2 + 19.*y1*y1 - 9.*y1*y2 + 9.*y2*y2)/(20.*Delta);
  Dek[2][3] = -27.*(x0*x0 - x0*x1 - x0*x2 + x1*x1 - x1*x2 + x2*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y1 - y1*y2 + y2*y2)/(20.*Delta);

  Dek[3][0] = -27.*(x0*x0 - x0*x1 - x0*x2 + x1*x1 - x1*x2 + x2*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y1 - y1*y2 + y2*y2)/(20.*Delta);
  Dek[3][1] = -27.*(x0*x0 - x0*x1 - x0*x2 + x1*x1 - x1*x2 + x2*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y1 - y1*y2 + y2*y2)/(20.*Delta);
  Dek[3][2] = -27.*(x0*x0 - x0*x1 - x0*x2 + x1*x1 - x1*x2 + x2*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y1 - y1*y2 + y2*y2)/(20.*Delta);
  Dek[3][3] =  81.*(x0*x0 - x0*x1 - x0*x2 + x1*x1 - x1*x2 + x2*x2 + y0*y0 - y0*y1 - y0*y2 + y1*y1 - y1*y2 + y2*y2)/(20.*Delta);

	return 0;
}

PetscErrorCode Calc3DP1_Dek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, PetscScalar **Dek)
{

	const PetscInt n0 = e_inter[en][0];
	const PetscInt n1 = e_inter[en][1];
	const PetscInt n2 = e_inter[en][2];
	const PetscInt n3 = e_inter[en][3];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
	const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];

	const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 + x1*y0*z2 - x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3 + x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0 + x3*y1*z2 + x3*y2*z0 - x3*y2*z1;

	Dek[0][0] = (x1*x1*y2*y2 - 2.*x1*x1*y2*y3 + x1*x1*y3*y3 + x1*x1*z2*z2 - 2.*x1*x1*z2*z3 + x1*x1*z3*z3 - 2.*x1*x2*y1*y2 + 2.*x1*x2*y1*y3 + 2.*x1*x2*y2*y3 - 2.*x1*x2*y3*y3 - 2.*x1*x2*z1*z2 + 2.*x1*x2*z1*z3 + 2.*x1*x2*z2*z3 - 2.*x1*x2*z3*z3 + 2.*x1*x3*y1*y2 - 2.*x1*x3*y1*y3 - 2.*x1*x3*y2*y2 + 2.*x1*x3*y2*y3 + 2.*x1*x3*z1*z2 - 2.*x1*x3*z1*z3 - 2.*x1*x3*z2*z2 + 2.*x1*x3*z2*z3 + x2*x2*y1*y1 - 2.*x2*x2*y1*y3 + x2*x2*y3*y3 + x2*x2*z1*z1 - 2.*x2*x2*z1*z3 + x2*x2*z3*z3 - 2.*x2*x3*y1*y1 + 2.*x2*x3*y1*y2 + 2.*x2*x3*y1*y3 - 2.*x2*x3*y2*y3 - 2.*x2*x3*z1*z1 + 2.*x2*x3*z1*z2 + 2.*x2*x3*z1*z3 - 2.*x2*x3*z2*z3 + x3*x3*y1*y1 - 2.*x3*x3*y1*y2 + x3*x3*y2*y2 + x3*x3*z1*z1 - 2.*x3*x3*z1*z2 + x3*x3*z2*z2 + y1*y1*z2*z2 - 2.*y1*y1*z2*z3 + y1*y1*z3*z3 - 2.*y1*y2*z1*z2 + 2.*y1*y2*z1*z3 + 2.*y1*y2*z2*z3 - 2.*y1*y2*z3*z3 + 2.*y1*y3*z1*z2 - 2.*y1*y3*z1*z3 - 2.*y1*y3*z2*z2 + 2.*y1*y3*z2*z3 + y2*y2*z1*z1 - 2.*y2*y2*z1*z3 + y2*y2*z3*z3 - 2.*y2*y3*z1*z1 + 2.*y2*y3*z1*z2 + 2.*y2*y3*z1*z3 - 2.*y2*y3*z2*z3 + y3*y3*z1*z1 - 2.*y3*y3*z1*z2 + y3*y3*z2*z2)/(6.*Delta);

	Dek[0][1] = -(x0*x1*y2*y2 - 2.*x0*x1*y2*y3 + x0*x1*y3*y3 + x0*x1*z2*z2 - 2.*x0*x1*z2*z3 + x0*x1*z3*z3 - x0*x2*y1*y2 + x0*x2*y1*y3 + x0*x2*y2*y3 - x0*x2*y3*y3 - x0*x2*z1*z2 + x0*x2*z1*z3 + x0*x2*z2*z3 - x0*x2*z3*z3 + x0*x3*y1*y2 - x0*x3*y1*y3 - x0*x3*y2*y2 + x0*x3*y2*y3 + x0*x3*z1*z2 - x0*x3*z1*z3 - x0*x3*z2*z2 + x0*x3*z2*z3 - x1*x2*y0*y2 + x1*x2*y0*y3 + x1*x2*y2*y3 - x1*x2*y3*y3 - x1*x2*z0*z2 + x1*x2*z0*z3 + x1*x2*z2*z3 - x1*x2*z3*z3 + x1*x3*y0*y2 - x1*x3*y0*y3 - x1*x3*y2*y2 + x1*x3*y2*y3 + x1*x3*z0*z2 - x1*x3*z0*z3 - x1*x3*z2*z2 + x1*x3*z2*z3 + x2*x2*y0*y1 - x2*x2*y0*y3 - x2*x2*y1*y3 + x2*x2*y3*y3 + x2*x2*z0*z1 - x2*x2*z0*z3 - x2*x2*z1*z3 + x2*x2*z3*z3 - 2.*x2*x3*y0*y1 + x2*x3*y0*y2 + x2*x3*y0*y3 + x2*x3*y1*y2 + x2*x3*y1*y3 - 2.*x2*x3*y2*y3 - 2.*x2*x3*z0*z1 + x2*x3*z0*z2 + x2*x3*z0*z3 + x2*x3*z1*z2 + x2*x3*z1*z3 - 2.*x2*x3*z2*z3 + x3*x3*y0*y1 - x3*x3*y0*y2 - x3*x3*y1*y2 + x3*x3*y2*y2 + x3*x3*z0*z1 - x3*x3*z0*z2 - x3*x3*z1*z2 + x3*x3*z2*z2 + y0*y1*z2*z2 - 2.*y0*y1*z2*z3 + y0*y1*z3*z3 - y0*y2*z1*z2 + y0*y2*z1*z3 + y0*y2*z2*z3 - y0*y2*z3*z3 + y0*y3*z1*z2 - y0*y3*z1*z3 - y0*y3*z2*z2 + y0*y3*z2*z3 - y1*y2*z0*z2 + y1*y2*z0*z3 + y1*y2*z2*z3 - y1*y2*z3*z3 + y1*y3*z0*z2 - y1*y3*z0*z3 - y1*y3*z2*z2 + y1*y3*z2*z3 + y2*y2*z0*z1 - y2*y2*z0*z3 - y2*y2*z1*z3 + y2*y2*z3*z3 - 2.*y2*y3*z0*z1 + y2*y3*z0*z2 + y2*y3*z0*z3 + y2*y3*z1*z2 + y2*y3*z1*z3 - 2.*y2*y3*z2*z3 + y3*y3*z0*z1 - y3*y3*z0*z2 - y3*y3*z1*z2 + y3*y3*z2*z2)/(6.*Delta);

	Dek[0][2] = (x0*x1*y1*y2 - x0*x1*y1*y3 - x0*x1*y2*y3 + x0*x1*y3*y3 + x0*x1*z1*z2 - x0*x1*z1*z3 - x0*x1*z2*z3 + x0*x1*z3*z3 - x0*x2*y1*y1 + 2.*x0*x2*y1*y3 - x0*x2*y3*y3 - x0*x2*z1*z1 + 2.*x0*x2*z1*z3 - x0*x2*z3*z3 + x0*x3*y1*y1 - x0*x3*y1*y2 - x0*x3*y1*y3 + x0*x3*y2*y3 + x0*x3*z1*z1 - x0*x3*z1*z2 - x0*x3*z1*z3 + x0*x3*z2*z3 - x1*x1*y0*y2 + x1*x1*y0*y3 + x1*x1*y2*y3 - x1*x1*y3*y3 - x1*x1*z0*z2 + x1*x1*z0*z3 + x1*x1*z2*z3 - x1*x1*z3*z3 + x1*x2*y0*y1 - x1*x2*y0*y3 - x1*x2*y1*y3 + x1*x2*y3*y3 + x1*x2*z0*z1 - x1*x2*z0*z3 - x1*x2*z1*z3 + x1*x2*z3*z3 - x1*x3*y0*y1 + 2.*x1*x3*y0*y2 - x1*x3*y0*y3 - x1*x3*y1*y2 + 2.*x1*x3*y1*y3 - x1*x3*y2*y3 - x1*x3*z0*z1 + 2.*x1*x3*z0*z2 - x1*x3*z0*z3 - x1*x3*z1*z2 + 2.*x1*x3*z1*z3 - x1*x3*z2*z3 - x2*x3*y0*y1 + x2*x3*y0*y3 + x2*x3*y1*y1 - x2*x3*y1*y3 - x2*x3*z0*z1 + x2*x3*z0*z3 + x2*x3*z1*z1 - x2*x3*z1*z3 + x3*x3*y0*y1 - x3*x3*y0*y2 - x3*x3*y1*y1 + x3*x3*y1*y2 + x3*x3*z0*z1 - x3*x3*z0*z2 - x3*x3*z1*z1 + x3*x3*z1*z2 + y0*y1*z1*z2 - y0*y1*z1*z3 - y0*y1*z2*z3 + y0*y1*z3*z3 - y0*y2*z1*z1 + 2.*y0*y2*z1*z3 - y0*y2*z3*z3 + y0*y3*z1*z1 - y0*y3*z1*z2 - y0*y3*z1*z3 + y0*y3*z2*z3 - y1*y1*z0*z2 + y1*y1*z0*z3 + y1*y1*z2*z3 - y1*y1*z3*z3 + y1*y2*z0*z1 - y1*y2*z0*z3 - y1*y2*z1*z3 + y1*y2*z3*z3 - y1*y3*z0*z1 + 2.*y1*y3*z0*z2 - y1*y3*z0*z3 - y1*y3*z1*z2 + 2.*y1*y3*z1*z3 - y1*y3*z2*z3 - y2*y3*z0*z1 + y2*y3*z0*z3 + y2*y3*z1*z1 - y2*y3*z1*z3 + y3*y3*z0*z1 - y3*y3*z0*z2 - y3*y3*z1*z1 + y3*y3*z1*z2)/(6.*Delta);

	Dek[0][3] = -(x0*x1*y1*y2 - x0*x1*y1*y3 - x0*x1*y2*y2 + x0*x1*y2*y3 + x0*x1*z1*z2 - x0*x1*z1*z3 - x0*x1*z2*z2 + x0*x1*z2*z3 - x0*x2*y1*y1 + x0*x2*y1*y2 + x0*x2*y1*y3 - x0*x2*y2*y3 - x0*x2*z1*z1 + x0*x2*z1*z2 + x0*x2*z1*z3 - x0*x2*z2*z3 + x0*x3*y1*y1 - 2.*x0*x3*y1*y2 + x0*x3*y2*y2 + x0*x3*z1*z1 - 2.*x0*x3*z1*z2 + x0*x3*z2*z2 - x1*x1*y0*y2 + x1*x1*y0*y3 + x1*x1*y2*y2 - x1*x1*y2*y3 - x1*x1*z0*z2 + x1*x1*z0*z3 + x1*x1*z2*z2 - x1*x1*z2*z3 + x1*x2*y0*y1 + x1*x2*y0*y2 - 2.*x1*x2*y0*y3 - 2.*x1*x2*y1*y2 + x1*x2*y1*y3 + x1*x2*y2*y3 + x1*x2*z0*z1 + x1*x2*z0*z2 - 2.*x1*x2*z0*z3 - 2.*x1*x2*z1*z2 + x1*x2*z1*z3 + x1*x2*z2*z3 - x1*x3*y0*y1 + x1*x3*y0*y2 + x1*x3*y1*y2 - x1*x3*y2*y2 - x1*x3*z0*z1 + x1*x3*z0*z2 + x1*x3*z1*z2 - x1*x3*z2*z2 - x2*x2*y0*y1 + x2*x2*y0*y3 + x2*x2*y1*y1 - x2*x2*y1*y3 - x2*x2*z0*z1 + x2*x2*z0*z3 + x2*x2*z1*z1 - x2*x2*z1*z3 + x2*x3*y0*y1 - x2*x3*y0*y2 - x2*x3*y1*y1 + x2*x3*y1*y2 + x2*x3*z0*z1 - x2*x3*z0*z2 - x2*x3*z1*z1 + x2*x3*z1*z2 + y0*y1*z1*z2 - y0*y1*z1*z3 - y0*y1*z2*z2 + y0*y1*z2*z3 - y0*y2*z1*z1 + y0*y2*z1*z2 + y0*y2*z1*z3 - y0*y2*z2*z3 + y0*y3*z1*z1 - 2.*y0*y3*z1*z2 + y0*y3*z2*z2 - y1*y1*z0*z2 + y1*y1*z0*z3 + y1*y1*z2*z2 - y1*y1*z2*z3 + y1*y2*z0*z1 + y1*y2*z0*z2 - 2.*y1*y2*z0*z3 - 2.*y1*y2*z1*z2 + y1*y2*z1*z3 + y1*y2*z2*z3 - y1*y3*z0*z1 + y1*y3*z0*z2 + y1*y3*z1*z2 - y1*y3*z2*z2 - y2*y2*z0*z1 + y2*y2*z0*z3 + y2*y2*z1*z1 - y2*y2*z1*z3 + y2*y3*z0*z1 - y2*y3*z0*z2 - y2*y3*z1*z1 + y2*y3*z1*z2)/(6.*Delta);

	Dek[1][0] = -(x0*x1*y2*y2 - 2.*x0*x1*y2*y3 + x0*x1*y3*y3 + x0*x1*z2*z2 - 2.*x0*x1*z2*z3 + x0*x1*z3*z3 - x0*x2*y1*y2 + x0*x2*y1*y3 + x0*x2*y2*y3 - x0*x2*y3*y3 - x0*x2*z1*z2 + x0*x2*z1*z3 + x0*x2*z2*z3 - x0*x2*z3*z3 + x0*x3*y1*y2 - x0*x3*y1*y3 - x0*x3*y2*y2 + x0*x3*y2*y3 + x0*x3*z1*z2 - x0*x3*z1*z3 - x0*x3*z2*z2 + x0*x3*z2*z3 - x1*x2*y0*y2 + x1*x2*y0*y3 + x1*x2*y2*y3 - x1*x2*y3*y3 - x1*x2*z0*z2 + x1*x2*z0*z3 + x1*x2*z2*z3 - x1*x2*z3*z3 + x1*x3*y0*y2 - x1*x3*y0*y3 - x1*x3*y2*y2 + x1*x3*y2*y3 + x1*x3*z0*z2 - x1*x3*z0*z3 - x1*x3*z2*z2 + x1*x3*z2*z3 + x2*x2*y0*y1 - x2*x2*y0*y3 - x2*x2*y1*y3 + x2*x2*y3*y3 + x2*x2*z0*z1 - x2*x2*z0*z3 - x2*x2*z1*z3 + x2*x2*z3*z3 - 2.*x2*x3*y0*y1 + x2*x3*y0*y2 + x2*x3*y0*y3 + x2*x3*y1*y2 + x2*x3*y1*y3 - 2.*x2*x3*y2*y3 - 2.*x2*x3*z0*z1 + x2*x3*z0*z2 + x2*x3*z0*z3 + x2*x3*z1*z2 + x2*x3*z1*z3 - 2.*x2*x3*z2*z3 + x3*x3*y0*y1 - x3*x3*y0*y2 - x3*x3*y1*y2 + x3*x3*y2*y2 + x3*x3*z0*z1 - x3*x3*z0*z2 - x3*x3*z1*z2 + x3*x3*z2*z2 + y0*y1*z2*z2 - 2.*y0*y1*z2*z3 + y0*y1*z3*z3 - y0*y2*z1*z2 + y0*y2*z1*z3 + y0*y2*z2*z3 - y0*y2*z3*z3 + y0*y3*z1*z2 - y0*y3*z1*z3 - y0*y3*z2*z2 + y0*y3*z2*z3 - y1*y2*z0*z2 + y1*y2*z0*z3 + y1*y2*z2*z3 - y1*y2*z3*z3 + y1*y3*z0*z2 - y1*y3*z0*z3 - y1*y3*z2*z2 + y1*y3*z2*z3 + y2*y2*z0*z1 - y2*y2*z0*z3 - y2*y2*z1*z3 + y2*y2*z3*z3 - 2.*y2*y3*z0*z1 + y2*y3*z0*z2 + y2*y3*z0*z3 + y2*y3*z1*z2 + y2*y3*z1*z3 - 2.*y2*y3*z2*z3 + y3*y3*z0*z1 - y3*y3*z0*z2 - y3*y3*z1*z2 + y3*y3*z2*z2)/(6.*Delta);

	Dek[1][1] = (x0*x0*y2*y2 - 2.*x0*x0*y2*y3 + x0*x0*y3*y3 + x0*x0*z2*z2 - 2.*x0*x0*z2*z3 + x0*x0*z3*z3 - 2.*x0*x2*y0*y2 + 2.*x0*x2*y0*y3 + 2.*x0*x2*y2*y3 - 2.*x0*x2*y3*y3 - 2.*x0*x2*z0*z2 + 2.*x0*x2*z0*z3 + 2.*x0*x2*z2*z3 - 2.*x0*x2*z3*z3 + 2.*x0*x3*y0*y2 - 2.*x0*x3*y0*y3 - 2.*x0*x3*y2*y2 + 2.*x0*x3*y2*y3 + 2.*x0*x3*z0*z2 - 2.*x0*x3*z0*z3 - 2.*x0*x3*z2*z2 + 2.*x0*x3*z2*z3 + x2*x2*y0*y0 - 2.*x2*x2*y0*y3 + x2*x2*y3*y3 + x2*x2*z0*z0 - 2.*x2*x2*z0*z3 + x2*x2*z3*z3 - 2.*x2*x3*y0*y0 + 2.*x2*x3*y0*y2 + 2.*x2*x3*y0*y3 - 2.*x2*x3*y2*y3 - 2.*x2*x3*z0*z0 + 2.*x2*x3*z0*z2 + 2.*x2*x3*z0*z3 - 2.*x2*x3*z2*z3 + x3*x3*y0*y0 - 2.*x3*x3*y0*y2 + x3*x3*y2*y2 + x3*x3*z0*z0 - 2.*x3*x3*z0*z2 + x3*x3*z2*z2 + y0*y0*z2*z2 - 2.*y0*y0*z2*z3 + y0*y0*z3*z3 - 2.*y0*y2*z0*z2 + 2.*y0*y2*z0*z3 + 2.*y0*y2*z2*z3 - 2.*y0*y2*z3*z3 + 2.*y0*y3*z0*z2 - 2.*y0*y3*z0*z3 - 2.*y0*y3*z2*z2 + 2.*y0*y3*z2*z3 + y2*y2*z0*z0 - 2.*y2*y2*z0*z3 + y2*y2*z3*z3 - 2.*y2*y3*z0*z0 + 2.*y2*y3*z0*z2 + 2.*y2*y3*z0*z3 - 2.*y2*y3*z2*z3 + y3*y3*z0*z0 - 2.*y3*y3*z0*z2 + y3*y3*z2*z2)/(6.*Delta);

	Dek[1][2] = -(x0*x0*y1*y2 - x0*x0*y1*y3 - x0*x0*y2*y3 + x0*x0*y3*y3 + x0*x0*z1*z2 - x0*x0*z1*z3 - x0*x0*z2*z3 + x0*x0*z3*z3 - x0*x1*y0*y2 + x0*x1*y0*y3 + x0*x1*y2*y3 - x0*x1*y3*y3 - x0*x1*z0*z2 + x0*x1*z0*z3 + x0*x1*z2*z3 - x0*x1*z3*z3 - x0*x2*y0*y1 + x0*x2*y0*y3 + x0*x2*y1*y3 - x0*x2*y3*y3 - x0*x2*z0*z1 + x0*x2*z0*z3 + x0*x2*z1*z3 - x0*x2*z3*z3 + x0*x3*y0*y1 + x0*x3*y0*y2 - 2.*x0*x3*y0*y3 - 2.*x0*x3*y1*y2 + x0*x3*y1*y3 + x0*x3*y2*y3 + x0*x3*z0*z1 + x0*x3*z0*z2 - 2.*x0*x3*z0*z3 - 2.*x0*x3*z1*z2 + x0*x3*z1*z3 + x0*x3*z2*z3 + x1*x2*y0*y0 - 2.*x1*x2*y0*y3 + x1*x2*y3*y3 + x1*x2*z0*z0 - 2.*x1*x2*z0*z3 + x1*x2*z3*z3 - x1*x3*y0*y0 + x1*x3*y0*y2 + x1*x3*y0*y3 - x1*x3*y2*y3 - x1*x3*z0*z0 + x1*x3*z0*z2 + x1*x3*z0*z3 - x1*x3*z2*z3 - x2*x3*y0*y0 + x2*x3*y0*y1 + x2*x3*y0*y3 - x2*x3*y1*y3 - x2*x3*z0*z0 + x2*x3*z0*z1 + x2*x3*z0*z3 - x2*x3*z1*z3 + x3*x3*y0*y0 - x3*x3*y0*y1 - x3*x3*y0*y2 + x3*x3*y1*y2 + x3*x3*z0*z0 - x3*x3*z0*z1 - x3*x3*z0*z2 + x3*x3*z1*z2 + y0*y0*z1*z2 - y0*y0*z1*z3 - y0*y0*z2*z3 + y0*y0*z3*z3 - y0*y1*z0*z2 + y0*y1*z0*z3 + y0*y1*z2*z3 - y0*y1*z3*z3 - y0*y2*z0*z1 + y0*y2*z0*z3 + y0*y2*z1*z3 - y0*y2*z3*z3 + y0*y3*z0*z1 + y0*y3*z0*z2 - 2.*y0*y3*z0*z3 - 2.*y0*y3*z1*z2 + y0*y3*z1*z3 + y0*y3*z2*z3 + y1*y2*z0*z0 - 2.*y1*y2*z0*z3 + y1*y2*z3*z3 - y1*y3*z0*z0 + y1*y3*z0*z2 + y1*y3*z0*z3 - y1*y3*z2*z3 - y2*y3*z0*z0 + y2*y3*z0*z1 + y2*y3*z0*z3 - y2*y3*z1*z3 + y3*y3*z0*z0 - y3*y3*z0*z1 - y3*y3*z0*z2 + y3*y3*z1*z2)/(6.*Delta);

	Dek[1][3] = (x0*x0*y1*y2 - x0*x0*y1*y3 - x0*x0*y2*y2 + x0*x0*y2*y3 + x0*x0*z1*z2 - x0*x0*z1*z3 - x0*x0*z2*z2 + x0*x0*z2*z3 - x0*x1*y0*y2 + x0*x1*y0*y3 + x0*x1*y2*y2 - x0*x1*y2*y3 - x0*x1*z0*z2 + x0*x1*z0*z3 + x0*x1*z2*z2 - x0*x1*z2*z3 - x0*x2*y0*y1 + 2.*x0*x2*y0*y2 - x0*x2*y0*y3 - x0*x2*y1*y2 + 2.*x0*x2*y1*y3 - x0*x2*y2*y3 - x0*x2*z0*z1 + 2.*x0*x2*z0*z2 - x0*x2*z0*z3 - x0*x2*z1*z2 + 2.*x0*x2*z1*z3 - x0*x2*z2*z3 + x0*x3*y0*y1 - x0*x3*y0*y2 - x0*x3*y1*y2 + x0*x3*y2*y2 + x0*x3*z0*z1 - x0*x3*z0*z2 - x0*x3*z1*z2 + x0*x3*z2*z2 + x1*x2*y0*y0 - x1*x2*y0*y2 - x1*x2*y0*y3 + x1*x2*y2*y3 + x1*x2*z0*z0 - x1*x2*z0*z2 - x1*x2*z0*z3 + x1*x2*z2*z3 - x1*x3*y0*y0 + 2.*x1*x3*y0*y2 - x1*x3*y2*y2 - x1*x3*z0*z0 + 2.*x1*x3*z0*z2 - x1*x3*z2*z2 - x2*x2*y0*y0 + x2*x2*y0*y1 + x2*x2*y0*y3 - x2*x2*y1*y3 - x2*x2*z0*z0 + x2*x2*z0*z1 + x2*x2*z0*z3 - x2*x2*z1*z3 + x2*x3*y0*y0 - x2*x3*y0*y1 - x2*x3*y0*y2 + x2*x3*y1*y2 + x2*x3*z0*z0 - x2*x3*z0*z1 - x2*x3*z0*z2 + x2*x3*z1*z2 + y0*y0*z1*z2 - y0*y0*z1*z3 - y0*y0*z2*z2 + y0*y0*z2*z3 - y0*y1*z0*z2 + y0*y1*z0*z3 + y0*y1*z2*z2 - y0*y1*z2*z3 - y0*y2*z0*z1 + 2.*y0*y2*z0*z2 - y0*y2*z0*z3 - y0*y2*z1*z2 + 2.*y0*y2*z1*z3 - y0*y2*z2*z3 + y0*y3*z0*z1 - y0*y3*z0*z2 - y0*y3*z1*z2 + y0*y3*z2*z2 + y1*y2*z0*z0 - y1*y2*z0*z2 - y1*y2*z0*z3 + y1*y2*z2*z3 - y1*y3*z0*z0 + 2.*y1*y3*z0*z2 - y1*y3*z2*z2 - y2*y2*z0*z0 + y2*y2*z0*z1 + y2*y2*z0*z3 - y2*y2*z1*z3 + y2*y3*z0*z0 - y2*y3*z0*z1 - y2*y3*z0*z2 + y2*y3*z1*z2)/(6.*Delta);

	Dek[2][0] = (x0*x1*y1*y2 - x0*x1*y1*y3 - x0*x1*y2*y3 + x0*x1*y3*y3 + x0*x1*z1*z2 - x0*x1*z1*z3 - x0*x1*z2*z3 + x0*x1*z3*z3 - x0*x2*y1*y1 + 2.*x0*x2*y1*y3 - x0*x2*y3*y3 - x0*x2*z1*z1 + 2.*x0*x2*z1*z3 - x0*x2*z3*z3 + x0*x3*y1*y1 - x0*x3*y1*y2 - x0*x3*y1*y3 + x0*x3*y2*y3 + x0*x3*z1*z1 - x0*x3*z1*z2 - x0*x3*z1*z3 + x0*x3*z2*z3 - x1*x1*y0*y2 + x1*x1*y0*y3 + x1*x1*y2*y3 - x1*x1*y3*y3 - x1*x1*z0*z2 + x1*x1*z0*z3 + x1*x1*z2*z3 - x1*x1*z3*z3 + x1*x2*y0*y1 - x1*x2*y0*y3 - x1*x2*y1*y3 + x1*x2*y3*y3 + x1*x2*z0*z1 - x1*x2*z0*z3 - x1*x2*z1*z3 + x1*x2*z3*z3 - x1*x3*y0*y1 + 2.*x1*x3*y0*y2 - x1*x3*y0*y3 - x1*x3*y1*y2 + 2.*x1*x3*y1*y3 - x1*x3*y2*y3 - x1*x3*z0*z1 + 2.*x1*x3*z0*z2 - x1*x3*z0*z3 - x1*x3*z1*z2 + 2.*x1*x3*z1*z3 - x1*x3*z2*z3 - x2*x3*y0*y1 + x2*x3*y0*y3 + x2*x3*y1*y1 - x2*x3*y1*y3 - x2*x3*z0*z1 + x2*x3*z0*z3 + x2*x3*z1*z1 - x2*x3*z1*z3 + x3*x3*y0*y1 - x3*x3*y0*y2 - x3*x3*y1*y1 + x3*x3*y1*y2 + x3*x3*z0*z1 - x3*x3*z0*z2 - x3*x3*z1*z1 + x3*x3*z1*z2 + y0*y1*z1*z2 - y0*y1*z1*z3 - y0*y1*z2*z3 + y0*y1*z3*z3 - y0*y2*z1*z1 + 2.*y0*y2*z1*z3 - y0*y2*z3*z3 + y0*y3*z1*z1 - y0*y3*z1*z2 - y0*y3*z1*z3 + y0*y3*z2*z3 - y1*y1*z0*z2 + y1*y1*z0*z3 + y1*y1*z2*z3 - y1*y1*z3*z3 + y1*y2*z0*z1 - y1*y2*z0*z3 - y1*y2*z1*z3 + y1*y2*z3*z3 - y1*y3*z0*z1 + 2.*y1*y3*z0*z2 - y1*y3*z0*z3 - y1*y3*z1*z2 + 2.*y1*y3*z1*z3 - y1*y3*z2*z3 - y2*y3*z0*z1 + y2*y3*z0*z3 + y2*y3*z1*z1 - y2*y3*z1*z3 + y3*y3*z0*z1 - y3*y3*z0*z2 - y3*y3*z1*z1 + y3*y3*z1*z2)/(6.*Delta);

	Dek[2][1] = -(x0*x0*y1*y2 - x0*x0*y1*y3 - x0*x0*y2*y3 + x0*x0*y3*y3 + x0*x0*z1*z2 - x0*x0*z1*z3 - x0*x0*z2*z3 + x0*x0*z3*z3 - x0*x1*y0*y2 + x0*x1*y0*y3 + x0*x1*y2*y3 - x0*x1*y3*y3 - x0*x1*z0*z2 + x0*x1*z0*z3 + x0*x1*z2*z3 - x0*x1*z3*z3 - x0*x2*y0*y1 + x0*x2*y0*y3 + x0*x2*y1*y3 - x0*x2*y3*y3 - x0*x2*z0*z1 + x0*x2*z0*z3 + x0*x2*z1*z3 - x0*x2*z3*z3 + x0*x3*y0*y1 + x0*x3*y0*y2 - 2.*x0*x3*y0*y3 - 2.*x0*x3*y1*y2 + x0*x3*y1*y3 + x0*x3*y2*y3 + x0*x3*z0*z1 + x0*x3*z0*z2 - 2.*x0*x3*z0*z3 - 2.*x0*x3*z1*z2 + x0*x3*z1*z3 + x0*x3*z2*z3 + x1*x2*y0*y0 - 2.*x1*x2*y0*y3 + x1*x2*y3*y3 + x1*x2*z0*z0 - 2.*x1*x2*z0*z3 + x1*x2*z3*z3 - x1*x3*y0*y0 + x1*x3*y0*y2 + x1*x3*y0*y3 - x1*x3*y2*y3 - x1*x3*z0*z0 + x1*x3*z0*z2 + x1*x3*z0*z3 - x1*x3*z2*z3 - x2*x3*y0*y0 + x2*x3*y0*y1 + x2*x3*y0*y3 - x2*x3*y1*y3 - x2*x3*z0*z0 + x2*x3*z0*z1 + x2*x3*z0*z3 - x2*x3*z1*z3 + x3*x3*y0*y0 - x3*x3*y0*y1 - x3*x3*y0*y2 + x3*x3*y1*y2 + x3*x3*z0*z0 - x3*x3*z0*z1 - x3*x3*z0*z2 + x3*x3*z1*z2 + y0*y0*z1*z2 - y0*y0*z1*z3 - y0*y0*z2*z3 + y0*y0*z3*z3 - y0*y1*z0*z2 + y0*y1*z0*z3 + y0*y1*z2*z3 - y0*y1*z3*z3 - y0*y2*z0*z1 + y0*y2*z0*z3 + y0*y2*z1*z3 - y0*y2*z3*z3 + y0*y3*z0*z1 + y0*y3*z0*z2 - 2.*y0*y3*z0*z3 - 2.*y0*y3*z1*z2 + y0*y3*z1*z3 + y0*y3*z2*z3 + y1*y2*z0*z0 - 2.*y1*y2*z0*z3 + y1*y2*z3*z3 - y1*y3*z0*z0 + y1*y3*z0*z2 + y1*y3*z0*z3 - y1*y3*z2*z3 - y2*y3*z0*z0 + y2*y3*z0*z1 + y2*y3*z0*z3 - y2*y3*z1*z3 + y3*y3*z0*z0 - y3*y3*z0*z1 - y3*y3*z0*z2 + y3*y3*z1*z2)/(6.*Delta);

	Dek[2][2] = (x0*x0*y1*y1 - 2.*x0*x0*y1*y3 + x0*x0*y3*y3 + x0*x0*z1*z1 - 2.*x0*x0*z1*z3 + x0*x0*z3*z3 - 2.*x0*x1*y0*y1 + 2.*x0*x1*y0*y3 + 2.*x0*x1*y1*y3 - 2.*x0*x1*y3*y3 - 2.*x0*x1*z0*z1 + 2.*x0*x1*z0*z3 + 2.*x0*x1*z1*z3 - 2.*x0*x1*z3*z3 + 2.*x0*x3*y0*y1 - 2.*x0*x3*y0*y3 - 2.*x0*x3*y1*y1 + 2.*x0*x3*y1*y3 + 2.*x0*x3*z0*z1 - 2.*x0*x3*z0*z3 - 2.*x0*x3*z1*z1 + 2.*x0*x3*z1*z3 + x1*x1*y0*y0 - 2.*x1*x1*y0*y3 + x1*x1*y3*y3 + x1*x1*z0*z0 - 2.*x1*x1*z0*z3 + x1*x1*z3*z3 - 2.*x1*x3*y0*y0 + 2.*x1*x3*y0*y1 + 2.*x1*x3*y0*y3 - 2.*x1*x3*y1*y3 - 2.*x1*x3*z0*z0 + 2.*x1*x3*z0*z1 + 2.*x1*x3*z0*z3 - 2.*x1*x3*z1*z3 + x3*x3*y0*y0 - 2.*x3*x3*y0*y1 + x3*x3*y1*y1 + x3*x3*z0*z0 - 2.*x3*x3*z0*z1 + x3*x3*z1*z1 + y0*y0*z1*z1 - 2.*y0*y0*z1*z3 + y0*y0*z3*z3 - 2.*y0*y1*z0*z1 + 2.*y0*y1*z0*z3 + 2.*y0*y1*z1*z3 - 2.*y0*y1*z3*z3 + 2.*y0*y3*z0*z1 - 2.*y0*y3*z0*z3 - 2.*y0*y3*z1*z1 + 2.*y0*y3*z1*z3 + y1*y1*z0*z0 - 2.*y1*y1*z0*z3 + y1*y1*z3*z3 - 2.*y1*y3*z0*z0 + 2.*y1*y3*z0*z1 + 2.*y1*y3*z0*z3 - 2.*y1*y3*z1*z3 + y3*y3*z0*z0 - 2.*y3*y3*z0*z1 + y3*y3*z1*z1)/(6.*Delta);

	Dek[2][3] = -(x0*x0*y1*y1 - x0*x0*y1*y2 - x0*x0*y1*y3 + x0*x0*y2*y3 + x0*x0*z1*z1 - x0*x0*z1*z2 - x0*x0*z1*z3 + x0*x0*z2*z3 - 2.*x0*x1*y0*y1 + x0*x1*y0*y2 + x0*x1*y0*y3 + x0*x1*y1*y2 + x0*x1*y1*y3 - 2.*x0*x1*y2*y3 - 2.*x0*x1*z0*z1 + x0*x1*z0*z2 + x0*x1*z0*z3 + x0*x1*z1*z2 + x0*x1*z1*z3 - 2.*x0*x1*z2*z3 + x0*x2*y0*y1 - x0*x2*y0*y3 - x0*x2*y1*y1 + x0*x2*y1*y3 + x0*x2*z0*z1 - x0*x2*z0*z3 - x0*x2*z1*z1 + x0*x2*z1*z3 + x0*x3*y0*y1 - x0*x3*y0*y2 - x0*x3*y1*y1 + x0*x3*y1*y2 + x0*x3*z0*z1 - x0*x3*z0*z2 - x0*x3*z1*z1 + x0*x3*z1*z2 + x1*x1*y0*y0 - x1*x1*y0*y2 - x1*x1*y0*y3 + x1*x1*y2*y3 + x1*x1*z0*z0 - x1*x1*z0*z2 - x1*x1*z0*z3 + x1*x1*z2*z3 - x1*x2*y0*y0 + x1*x2*y0*y1 + x1*x2*y0*y3 - x1*x2*y1*y3 - x1*x2*z0*z0 + x1*x2*z0*z1 + x1*x2*z0*z3 - x1*x2*z1*z3 - x1*x3*y0*y0 + x1*x3*y0*y1 + x1*x3*y0*y2 - x1*x3*y1*y2 - x1*x3*z0*z0 + x1*x3*z0*z1 + x1*x3*z0*z2 - x1*x3*z1*z2 + x2*x3*y0*y0 - 2.*x2*x3*y0*y1 + x2*x3*y1*y1 + x2*x3*z0*z0 - 2.*x2*x3*z0*z1 + x2*x3*z1*z1 + y0*y0*z1*z1 - y0*y0*z1*z2 - y0*y0*z1*z3 + y0*y0*z2*z3 - 2.*y0*y1*z0*z1 + y0*y1*z0*z2 + y0*y1*z0*z3 + y0*y1*z1*z2 + y0*y1*z1*z3 - 2.*y0*y1*z2*z3 + y0*y2*z0*z1 - y0*y2*z0*z3 - y0*y2*z1*z1 + y0*y2*z1*z3 + y0*y3*z0*z1 - y0*y3*z0*z2 - y0*y3*z1*z1 + y0*y3*z1*z2 + y1*y1*z0*z0 - y1*y1*z0*z2 - y1*y1*z0*z3 + y1*y1*z2*z3 - y1*y2*z0*z0 + y1*y2*z0*z1 + y1*y2*z0*z3 - y1*y2*z1*z3 - y1*y3*z0*z0 + y1*y3*z0*z1 + y1*y3*z0*z2 - y1*y3*z1*z2 + y2*y3*z0*z0 - 2.*y2*y3*z0*z1 + y2*y3*z1*z1)/(6.*Delta);

	Dek[3][0] = -(x0*x1*y1*y2 - x0*x1*y1*y3 - x0*x1*y2*y2 + x0*x1*y2*y3 + x0*x1*z1*z2 - x0*x1*z1*z3 - x0*x1*z2*z2 + x0*x1*z2*z3 - x0*x2*y1*y1 + x0*x2*y1*y2 + x0*x2*y1*y3 - x0*x2*y2*y3 - x0*x2*z1*z1 + x0*x2*z1*z2 + x0*x2*z1*z3 - x0*x2*z2*z3 + x0*x3*y1*y1 - 2.*x0*x3*y1*y2 + x0*x3*y2*y2 + x0*x3*z1*z1 - 2.*x0*x3*z1*z2 + x0*x3*z2*z2 - x1*x1*y0*y2 + x1*x1*y0*y3 + x1*x1*y2*y2 - x1*x1*y2*y3 - x1*x1*z0*z2 + x1*x1*z0*z3 + x1*x1*z2*z2 - x1*x1*z2*z3 + x1*x2*y0*y1 + x1*x2*y0*y2 - 2.*x1*x2*y0*y3 - 2.*x1*x2*y1*y2 + x1*x2*y1*y3 + x1*x2*y2*y3 + x1*x2*z0*z1 + x1*x2*z0*z2 - 2.*x1*x2*z0*z3 - 2.*x1*x2*z1*z2 + x1*x2*z1*z3 + x1*x2*z2*z3 - x1*x3*y0*y1 + x1*x3*y0*y2 + x1*x3*y1*y2 - x1*x3*y2*y2 - x1*x3*z0*z1 + x1*x3*z0*z2 + x1*x3*z1*z2 - x1*x3*z2*z2 - x2*x2*y0*y1 + x2*x2*y0*y3 + x2*x2*y1*y1 - x2*x2*y1*y3 - x2*x2*z0*z1 + x2*x2*z0*z3 + x2*x2*z1*z1 - x2*x2*z1*z3 + x2*x3*y0*y1 - x2*x3*y0*y2 - x2*x3*y1*y1 + x2*x3*y1*y2 + x2*x3*z0*z1 - x2*x3*z0*z2 - x2*x3*z1*z1 + x2*x3*z1*z2 + y0*y1*z1*z2 - y0*y1*z1*z3 - y0*y1*z2*z2 + y0*y1*z2*z3 - y0*y2*z1*z1 + y0*y2*z1*z2 + y0*y2*z1*z3 - y0*y2*z2*z3 + y0*y3*z1*z1 - 2.*y0*y3*z1*z2 + y0*y3*z2*z2 - y1*y1*z0*z2 + y1*y1*z0*z3 + y1*y1*z2*z2 - y1*y1*z2*z3 + y1*y2*z0*z1 + y1*y2*z0*z2 - 2.*y1*y2*z0*z3 - 2.*y1*y2*z1*z2 + y1*y2*z1*z3 + y1*y2*z2*z3 - y1*y3*z0*z1 + y1*y3*z0*z2 + y1*y3*z1*z2 - y1*y3*z2*z2 - y2*y2*z0*z1 + y2*y2*z0*z3 + y2*y2*z1*z1 - y2*y2*z1*z3 + y2*y3*z0*z1 - y2*y3*z0*z2 - y2*y3*z1*z1 + y2*y3*z1*z2)/(6.*Delta);

	Dek[3][1] = (x0*x0*y1*y2 - x0*x0*y1*y3 - x0*x0*y2*y2 + x0*x0*y2*y3 + x0*x0*z1*z2 - x0*x0*z1*z3 - x0*x0*z2*z2 + x0*x0*z2*z3 - x0*x1*y0*y2 + x0*x1*y0*y3 + x0*x1*y2*y2 - x0*x1*y2*y3 - x0*x1*z0*z2 + x0*x1*z0*z3 + x0*x1*z2*z2 - x0*x1*z2*z3 - x0*x2*y0*y1 + 2.*x0*x2*y0*y2 - x0*x2*y0*y3 - x0*x2*y1*y2 + 2.*x0*x2*y1*y3 - x0*x2*y2*y3 - x0*x2*z0*z1 + 2.*x0*x2*z0*z2 - x0*x2*z0*z3 - x0*x2*z1*z2 + 2.*x0*x2*z1*z3 - x0*x2*z2*z3 + x0*x3*y0*y1 - x0*x3*y0*y2 - x0*x3*y1*y2 + x0*x3*y2*y2 + x0*x3*z0*z1 - x0*x3*z0*z2 - x0*x3*z1*z2 + x0*x3*z2*z2 + x1*x2*y0*y0 - x1*x2*y0*y2 - x1*x2*y0*y3 + x1*x2*y2*y3 + x1*x2*z0*z0 - x1*x2*z0*z2 - x1*x2*z0*z3 + x1*x2*z2*z3 - x1*x3*y0*y0 + 2.*x1*x3*y0*y2 - x1*x3*y2*y2 - x1*x3*z0*z0 + 2.*x1*x3*z0*z2 - x1*x3*z2*z2 - x2*x2*y0*y0 + x2*x2*y0*y1 + x2*x2*y0*y3 - x2*x2*y1*y3 - x2*x2*z0*z0 + x2*x2*z0*z1 + x2*x2*z0*z3 - x2*x2*z1*z3 + x2*x3*y0*y0 - x2*x3*y0*y1 - x2*x3*y0*y2 + x2*x3*y1*y2 + x2*x3*z0*z0 - x2*x3*z0*z1 - x2*x3*z0*z2 + x2*x3*z1*z2 + y0*y0*z1*z2 - y0*y0*z1*z3 - y0*y0*z2*z2 + y0*y0*z2*z3 - y0*y1*z0*z2 + y0*y1*z0*z3 + y0*y1*z2*z2 - y0*y1*z2*z3 - y0*y2*z0*z1 + 2.*y0*y2*z0*z2 - y0*y2*z0*z3 - y0*y2*z1*z2 + 2.*y0*y2*z1*z3 - y0*y2*z2*z3 + y0*y3*z0*z1 - y0*y3*z0*z2 - y0*y3*z1*z2 + y0*y3*z2*z2 + y1*y2*z0*z0 - y1*y2*z0*z2 - y1*y2*z0*z3 + y1*y2*z2*z3 - y1*y3*z0*z0 + 2.*y1*y3*z0*z2 - y1*y3*z2*z2 - y2*y2*z0*z0 + y2*y2*z0*z1 + y2*y2*z0*z3 - y2*y2*z1*z3 + y2*y3*z0*z0 - y2*y3*z0*z1 - y2*y3*z0*z2 + y2*y3*z1*z2)/(6.*Delta);

	Dek[3][2] = -(x0*x0*y1*y1 - x0*x0*y1*y2 - x0*x0*y1*y3 + x0*x0*y2*y3 + x0*x0*z1*z1 - x0*x0*z1*z2 - x0*x0*z1*z3 + x0*x0*z2*z3 - 2.*x0*x1*y0*y1 + x0*x1*y0*y2 + x0*x1*y0*y3 + x0*x1*y1*y2 + x0*x1*y1*y3 - 2.*x0*x1*y2*y3 - 2.*x0*x1*z0*z1 + x0*x1*z0*z2 + x0*x1*z0*z3 + x0*x1*z1*z2 + x0*x1*z1*z3 - 2.*x0*x1*z2*z3 + x0*x2*y0*y1 - x0*x2*y0*y3 - x0*x2*y1*y1 + x0*x2*y1*y3 + x0*x2*z0*z1 - x0*x2*z0*z3 - x0*x2*z1*z1 + x0*x2*z1*z3 + x0*x3*y0*y1 - x0*x3*y0*y2 - x0*x3*y1*y1 + x0*x3*y1*y2 + x0*x3*z0*z1 - x0*x3*z0*z2 - x0*x3*z1*z1 + x0*x3*z1*z2 + x1*x1*y0*y0 - x1*x1*y0*y2 - x1*x1*y0*y3 + x1*x1*y2*y3 + x1*x1*z0*z0 - x1*x1*z0*z2 - x1*x1*z0*z3 + x1*x1*z2*z3 - x1*x2*y0*y0 + x1*x2*y0*y1 + x1*x2*y0*y3 - x1*x2*y1*y3 - x1*x2*z0*z0 + x1*x2*z0*z1 + x1*x2*z0*z3 - x1*x2*z1*z3 - x1*x3*y0*y0 + x1*x3*y0*y1 + x1*x3*y0*y2 - x1*x3*y1*y2 - x1*x3*z0*z0 + x1*x3*z0*z1 + x1*x3*z0*z2 - x1*x3*z1*z2 + x2*x3*y0*y0 - 2.*x2*x3*y0*y1 + x2*x3*y1*y1 + x2*x3*z0*z0 - 2.*x2*x3*z0*z1 + x2*x3*z1*z1 + y0*y0*z1*z1 - y0*y0*z1*z2 - y0*y0*z1*z3 + y0*y0*z2*z3 - 2.*y0*y1*z0*z1 + y0*y1*z0*z2 + y0*y1*z0*z3 + y0*y1*z1*z2 + y0*y1*z1*z3 - 2.*y0*y1*z2*z3 + y0*y2*z0*z1 - y0*y2*z0*z3 - y0*y2*z1*z1 + y0*y2*z1*z3 + y0*y3*z0*z1 - y0*y3*z0*z2 - y0*y3*z1*z1 + y0*y3*z1*z2 + y1*y1*z0*z0 - y1*y1*z0*z2 - y1*y1*z0*z3 + y1*y1*z2*z3 - y1*y2*z0*z0 + y1*y2*z0*z1 + y1*y2*z0*z3 - y1*y2*z1*z3 - y1*y3*z0*z0 + y1*y3*z0*z1 + y1*y3*z0*z2 - y1*y3*z1*z2 + y2*y3*z0*z0 - 2.*y2*y3*z0*z1 + y2*y3*z1*z1)/(6.*Delta);

	Dek[3][3] = (x0*x0*y1*y1 - 2.*x0*x0*y1*y2 + x0*x0*y2*y2 + x0*x0*z1*z1 - 2.*x0*x0*z1*z2 + x0*x0*z2*z2 - 2.*x0*x1*y0*y1 + 2.*x0*x1*y0*y2 + 2.*x0*x1*y1*y2 - 2.*x0*x1*y2*y2 - 2.*x0*x1*z0*z1 + 2.*x0*x1*z0*z2 + 2.*x0*x1*z1*z2 - 2.*x0*x1*z2*z2 + 2.*x0*x2*y0*y1 - 2.*x0*x2*y0*y2 - 2.*x0*x2*y1*y1 + 2.*x0*x2*y1*y2 + 2.*x0*x2*z0*z1 - 2.*x0*x2*z0*z2 - 2.*x0*x2*z1*z1 + 2.*x0*x2*z1*z2 + x1*x1*y0*y0 - 2.*x1*x1*y0*y2 + x1*x1*y2*y2 + x1*x1*z0*z0 - 2.*x1*x1*z0*z2 + x1*x1*z2*z2 - 2.*x1*x2*y0*y0 + 2.*x1*x2*y0*y1 + 2.*x1*x2*y0*y2 - 2.*x1*x2*y1*y2 - 2.*x1*x2*z0*z0 + 2.*x1*x2*z0*z1 + 2.*x1*x2*z0*z2 - 2.*x1*x2*z1*z2 + x2*x2*y0*y0 - 2.*x2*x2*y0*y1 + x2*x2*y1*y1 + x2*x2*z0*z0 - 2.*x2*x2*z0*z1 + x2*x2*z1*z1 + y0*y0*z1*z1 - 2.*y0*y0*z1*z2 + y0*y0*z2*z2 - 2.*y0*y1*z0*z1 + 2.*y0*y1*z0*z2 + 2.*y0*y1*z1*z2 - 2.*y0*y1*z2*z2 + 2.*y0*y2*z0*z1 - 2.*y0*y2*z0*z2 - 2.*y0*y2*z1*z1 + 2.*y0*y2*z1*z2 + y1*y1*z0*z0 - 2.*y1*y1*z0*z2 + y1*y1*z2*z2 - 2.*y1*y2*z0*z0 + 2.*y1*y2*z0*z1 + 2.*y1*y2*z0*z2 - 2.*y1*y2*z1*z2 + y2*y2*z0*z0 - 2.*y2*y2*z0*z1 + y2*y2*z1*z1)/(6.*Delta);

	return 0;
}


PetscErrorCode Calc3DP2_Dek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **eP2_inter, PetscScalar **Dek)
{
	const PetscInt n0 = eP2_inter[en][0];
	const PetscInt n1 = eP2_inter[en][1];
	const PetscInt n2 = eP2_inter[en][2];
	const PetscInt n3 = eP2_inter[en][3];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
	const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];
	
	const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 + x1*y0*z2 
            -x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3  
            +x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0  
            +x3*y1*z2 + x3*y2*z0 - x3*y2*z1;

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

	Dek[0][0] =  (C0xy*C0xy + C0xz*C0xz + C0yz*C0yz)/(10.*Delta);
	Dek[0][1] = -(C0xy*C1xy + C0xz*C1xz + C0yz*C1yz)/(30.*Delta);
	Dek[0][2] = -(C0xy*C2xy + C0xz*C2xz + C0yz*C2yz)/(30.*Delta);
	Dek[0][3] = -(C0xy*C3xy + C0xz*C3xz + C0yz*C3yz)/(30.*Delta);
	Dek[0][4] = -(C0xy*C0xy - 3.*C0xy*C1xy + C0xz*C0xz - 3.*C0xz*C1xz + C0yz*C0yz - 3.*C0yz*C1yz)/(30.*Delta);
	Dek[0][5] = -(C0xy*C1xy + C0xy*C2xy + C0xz*C1xz + C0xz*C2xz + C0yz*C1yz + C0yz*C2yz)/(30.*Delta);
	Dek[0][6] = -(C0xy*C0xy - 3.*C0xy*C2xy + C0xz*C0xz - 3.*C0xz*C2xz + C0yz*C0yz - 3.*C0yz*C2yz)/(30.*Delta);
	Dek[0][7] = -(C0xy*C0xy - 3.*C0xy*C3xy + C0xz*C0xz - 3.*C0xz*C3xz + C0yz*C0yz - 3.*C0yz*C3yz)/(30.*Delta);
	Dek[0][8] = -(C0xy*C1xy + C0xy*C3xy + C0xz*C1xz + C0xz*C3xz + C0yz*C1yz + C0yz*C3yz)/(30.*Delta);
	Dek[0][9] = -(C0xy*C2xy + C0xy*C3xy + C0xz*C2xz + C0xz*C3xz + C0yz*C2yz + C0yz*C3yz)/(30.*Delta);
	Dek[1][0] = -(C0xy*C1xy + C0xz*C1xz + C0yz*C1yz)/(30.*Delta);
	Dek[1][1] =  (C1xy*C1xy + C1xz*C1xz + C1yz*C1yz)/(10.*Delta);
	Dek[1][2] = -(C1xy*C2xy + C1xz*C2xz + C1yz*C2yz)/(30.*Delta);
	Dek[1][3] = -(C1xy*C3xy + C1xz*C3xz + C1yz*C3yz)/(30.*Delta);
	Dek[1][4] = (3*C0xy*C1xy + 3.*C0xz*C1xz + 3.*C0yz*C1yz - C1xy*C1xy - C1xz*C1xz - C1yz*C1yz)/(30.*Delta);
	Dek[1][5] = -(C1xy*C1xy - 3.*C1xy*C2xy + C1xz*C1xz - 3.*C1xz*C2xz + C1yz*C1yz - 3.*C1yz*C2yz)/(30.*Delta);
	Dek[1][6] = -(C0xy*C1xy + C0xz*C1xz + C0yz*C1yz + C1xy*C2xy + C1xz*C2xz + C1yz*C2yz)/(30.*Delta);
	Dek[1][7] = -(C0xy*C1xy + C0xz*C1xz + C0yz*C1yz + C1xy*C3xy + C1xz*C3xz + C1yz*C3yz)/(30.*Delta);
	Dek[1][8] = -(C1xy*C1xy - 3.*C1xy*C3xy + C1xz*C1xz - 3.*C1xz*C3xz + C1yz*C1yz - 3.*C1yz*C3yz)/(30.*Delta);
	Dek[1][9] = -(C1xy*C2xy + C1xy*C3xy + C1xz*C2xz + C1xz*C3xz + C1yz*C2yz + C1yz*C3yz)/(30.*Delta);
	Dek[2][0] = -(C0xy*C2xy + C0xz*C2xz + C0yz*C2yz)/(30.*Delta);
	Dek[2][1] = -(C1xy*C2xy + C1xz*C2xz + C1yz*C2yz)/(30.*Delta);
	Dek[2][2] = (C2xy*C2xy + C2xz*C2xz + C2yz*C2yz)/(10.*Delta);
	Dek[2][3] = -(C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(30.*Delta);
	Dek[2][4] = -(C0xy*C2xy + C0xz*C2xz + C0yz*C2yz + C1xy*C2xy + C1xz*C2xz + C1yz*C2yz)/(30.*Delta);
	Dek[2][5] = (3*C1xy*C2xy + 3.*C1xz*C2xz + 3.*C1yz*C2yz - C2xy*C2xy - C2xz*C2xz - C2yz*C2yz)/(30.*Delta);
	Dek[2][6] = (3*C0xy*C2xy + 3.*C0xz*C2xz + 3.*C0yz*C2yz - C2xy*C2xy - C2xz*C2xz - C2yz*C2yz)/(30.*Delta);
	Dek[2][7] = -(C0xy*C2xy + C0xz*C2xz + C0yz*C2yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(30.*Delta);
	Dek[2][8] = -(C1xy*C2xy + C1xz*C2xz + C1yz*C2yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(30.*Delta);
	Dek[2][9] = -(C2xy*C2xy - 3.*C2xy*C3xy + C2xz*C2xz - 3.*C2xz*C3xz + C2yz*C2yz - 3.*C2yz*C3yz)/(30.*Delta);
	Dek[3][0] = -(C0xy*C3xy + C0xz*C3xz + C0yz*C3yz)/(30.*Delta);
	Dek[3][1] = -(C1xy*C3xy + C1xz*C3xz + C1yz*C3yz)/(30.*Delta);
	Dek[3][2] = -(C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(30.*Delta);
	Dek[3][3] = (C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(10.*Delta);
	Dek[3][4] = -(C0xy*C3xy + C0xz*C3xz + C0yz*C3yz + C1xy*C3xy + C1xz*C3xz + C1yz*C3yz)/(30.*Delta);
	Dek[3][5] = -(C1xy*C3xy + C1xz*C3xz + C1yz*C3yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(30.*Delta);
	Dek[3][6] = -(C0xy*C3xy + C0xz*C3xz + C0yz*C3yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(30.*Delta);
	Dek[3][7] = (3*C0xy*C3xy + 3.*C0xz*C3xz + 3.*C0yz*C3yz - C3xy*C3xy - C3xz*C3xz - C3yz*C3yz)/(30.*Delta);
	Dek[3][8] = (3*C1xy*C3xy + 3.*C1xz*C3xz + 3.*C1yz*C3yz - C3xy*C3xy - C3xz*C3xz - C3yz*C3yz)/(30.*Delta);
	Dek[3][9] = (3*C2xy*C3xy + 3.*C2xz*C3xz + 3.*C2yz*C3yz - C3xy*C3xy - C3xz*C3xz - C3yz*C3yz)/(30.*Delta);
	Dek[4][0] = -(C0xy*C0xy - 3.*C0xy*C1xy + C0xz*C0xz - 3.*C0xz*C1xz + C0yz*C0yz - 3.*C0yz*C1yz)/(30.*Delta);
	Dek[4][1] = (3*C0xy*C1xy + 3.*C0xz*C1xz + 3.*C0yz*C1yz - C1xy*C1xy - C1xz*C1xz - C1yz*C1yz)/(30.*Delta);
	Dek[4][2] = -(C0xy*C2xy + C0xz*C2xz + C0yz*C2yz + C1xy*C2xy + C1xz*C2xz + C1yz*C2yz)/(30.*Delta);
	Dek[4][3] = -(C0xy*C3xy + C0xz*C3xz + C0yz*C3yz + C1xy*C3xy + C1xz*C3xz + C1yz*C3yz)/(30.*Delta);
	Dek[4][4] = 4.*(C0xy*C0xy + C0xy*C1xy + C0xz*C0xz + C0xz*C1xz + C0yz*C0yz + C0yz*C1yz + C1xy*C1xy + C1xz*C1xz + C1yz*C1yz)/(15.*Delta);
	Dek[4][5] = 2.*(C0xy*C1xy + 2.*C0xy*C2xy + C0xz*C1xz + 2.*C0xz*C2xz + C0yz*C1yz + 2.*C0yz*C2yz + C1xy*C1xy + C1xy*C2xy + C1xz*C1xz + C1xz*C2xz + C1yz*C1yz + C1yz*C2yz)/(15.*Delta);
	Dek[4][6] = 2.*(C0xy*C0xy + C0xy*C1xy + C0xy*C2xy + C0xz*C0xz + C0xz*C1xz + C0xz*C2xz + C0yz*C0yz + C0yz*C1yz + C0yz*C2yz + 2.*C1xy*C2xy + 2.*C1xz*C2xz + 2.*C1yz*C2yz)/(15.*Delta);
	Dek[4][7] = 2.*(C0xy*C0xy + C0xy*C1xy + C0xy*C3xy + C0xz*C0xz + C0xz*C1xz + C0xz*C3xz + C0yz*C0yz + C0yz*C1yz + C0yz*C3yz + 2.*C1xy*C3xy + 2.*C1xz*C3xz + 2.*C1yz*C3yz)/(15.*Delta);
	Dek[4][8] = 2.*(C0xy*C1xy + 2.*C0xy*C3xy + C0xz*C1xz + 2.*C0xz*C3xz + C0yz*C1yz + 2.*C0yz*C3yz + C1xy*C1xy + C1xy*C3xy + C1xz*C1xz + C1xz*C3xz + C1yz*C1yz + C1yz*C3yz)/(15.*Delta);
	Dek[4][9] = 2.*(C0xy*C2xy + C0xy*C3xy + C0xz*C2xz + C0xz*C3xz + C0yz*C2yz + C0yz*C3yz + C1xy*C2xy + C1xy*C3xy + C1xz*C2xz + C1xz*C3xz + C1yz*C2yz + C1yz*C3yz)/(15.*Delta);
	Dek[5][0] = -(C0xy*C1xy + C0xy*C2xy + C0xz*C1xz + C0xz*C2xz + C0yz*C1yz + C0yz*C2yz)/(30.*Delta);
	Dek[5][1] = -(C1xy*C1xy - 3.*C1xy*C2xy + C1xz*C1xz - 3.*C1xz*C2xz + C1yz*C1yz - 3.*C1yz*C2yz)/(30.*Delta);
	Dek[5][2] = (3*C1xy*C2xy + 3.*C1xz*C2xz + 3.*C1yz*C2yz - C2xy*C2xy - C2xz*C2xz - C2yz*C2yz)/(30.*Delta);
	Dek[5][3] = -(C1xy*C3xy + C1xz*C3xz + C1yz*C3yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(30.*Delta);
	Dek[5][4] = 2.*(C0xy*C1xy + 2.*C0xy*C2xy + C0xz*C1xz + 2.*C0xz*C2xz + C0yz*C1yz + 2.*C0yz*C2yz + C1xy*C1xy + C1xy*C2xy + C1xz*C1xz + C1xz*C2xz + C1yz*C1yz + C1yz*C2yz)/(15.*Delta);
	Dek[5][5] = 4.*(C1xy*C1xy + C1xy*C2xy + C1xz*C1xz + C1xz*C2xz + C1yz*C1yz + C1yz*C2yz + C2xy*C2xy + C2xz*C2xz + C2yz*C2yz)/(15.*Delta);
	Dek[5][6] = 2.*(2*C0xy*C1xy + C0xy*C2xy + 2.*C0xz*C1xz + C0xz*C2xz + 2.*C0yz*C1yz + C0yz*C2yz + C1xy*C2xy + C1xz*C2xz + C1yz*C2yz + C2xy*C2xy + C2xz*C2xz + C2yz*C2yz)/(15.*Delta);
	Dek[5][7] = 2.*(C0xy*C1xy + C0xy*C2xy + C0xz*C1xz + C0xz*C2xz + C0yz*C1yz + C0yz*C2yz + C1xy*C3xy + C1xz*C3xz + C1yz*C3yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(15.*Delta);
	Dek[5][8] = 2.*(C1xy*C1xy + C1xy*C2xy + C1xy*C3xy + C1xz*C1xz + C1xz*C2xz + C1xz*C3xz + C1yz*C1yz + C1yz*C2yz + C1yz*C3yz + 2.*C2xy*C3xy + 2.*C2xz*C3xz + 2.*C2yz*C3yz)/(15.*Delta);
	Dek[5][9] = 2.*(C1xy*C2xy + 2.*C1xy*C3xy + C1xz*C2xz + 2.*C1xz*C3xz + C1yz*C2yz + 2.*C1yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz)/(15.*Delta);
	Dek[6][0] = -(C0xy*C0xy - 3.*C0xy*C2xy + C0xz*C0xz - 3.*C0xz*C2xz + C0yz*C0yz - 3.*C0yz*C2yz)/(30.*Delta);
	Dek[6][1] = -(C0xy*C1xy + C0xz*C1xz + C0yz*C1yz + C1xy*C2xy + C1xz*C2xz + C1yz*C2yz)/(30.*Delta);
	Dek[6][2] = (3*C0xy*C2xy + 3.*C0xz*C2xz + 3.*C0yz*C2yz - C2xy*C2xy - C2xz*C2xz - C2yz*C2yz)/(30.*Delta);
	Dek[6][3] = -(C0xy*C3xy + C0xz*C3xz + C0yz*C3yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(30.*Delta);
	Dek[6][4] = 2.*(C0xy*C0xy + C0xy*C1xy + C0xy*C2xy + C0xz*C0xz + C0xz*C1xz + C0xz*C2xz + C0yz*C0yz + C0yz*C1yz + C0yz*C2yz + 2.*C1xy*C2xy + 2.*C1xz*C2xz + 2.*C1yz*C2yz)/(15.*Delta);
	Dek[6][5] = 2.*(2*C0xy*C1xy + C0xy*C2xy + 2.*C0xz*C1xz + C0xz*C2xz + 2.*C0yz*C1yz + C0yz*C2yz + C1xy*C2xy + C1xz*C2xz + C1yz*C2yz + C2xy*C2xy + C2xz*C2xz + C2yz*C2yz)/(15.*Delta);
	Dek[6][6] = 4.*(C0xy*C0xy + C0xy*C2xy + C0xz*C0xz + C0xz*C2xz + C0yz*C0yz + C0yz*C2yz + C2xy*C2xy + C2xz*C2xz + C2yz*C2yz)/(15.*Delta);
	Dek[6][7] = 2.*(C0xy*C0xy + C0xy*C2xy + C0xy*C3xy + C0xz*C0xz + C0xz*C2xz + C0xz*C3xz + C0yz*C0yz + C0yz*C2yz + C0yz*C3yz + 2.*C2xy*C3xy + 2.*C2xz*C3xz + 2.*C2yz*C3yz)/(15.*Delta);
	Dek[6][8] = 2.*(C0xy*C1xy + C0xy*C3xy + C0xz*C1xz + C0xz*C3xz + C0yz*C1yz + C0yz*C3yz + C1xy*C2xy + C1xz*C2xz + C1yz*C2yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(15.*Delta);
	Dek[6][9] = 2.*(C0xy*C2xy + 2.*C0xy*C3xy + C0xz*C2xz + 2.*C0xz*C3xz + C0yz*C2yz + 2.*C0yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz)/(15.*Delta);
	Dek[7][0] = -(C0xy*C0xy - 3.*C0xy*C3xy + C0xz*C0xz - 3.*C0xz*C3xz + C0yz*C0yz - 3.*C0yz*C3yz)/(30.*Delta);
	Dek[7][1] = -(C0xy*C1xy + C0xz*C1xz + C0yz*C1yz + C1xy*C3xy + C1xz*C3xz + C1yz*C3yz)/(30.*Delta);
	Dek[7][2] = -(C0xy*C2xy + C0xz*C2xz + C0yz*C2yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(30.*Delta);
	Dek[7][3] = (3*C0xy*C3xy + 3.*C0xz*C3xz + 3.*C0yz*C3yz - C3xy*C3xy - C3xz*C3xz - C3yz*C3yz)/(30.*Delta);
	Dek[7][4] = 2.*(C0xy*C0xy + C0xy*C1xy + C0xy*C3xy + C0xz*C0xz + C0xz*C1xz + C0xz*C3xz + C0yz*C0yz + C0yz*C1yz + C0yz*C3yz + 2.*C1xy*C3xy + 2.*C1xz*C3xz + 2.*C1yz*C3yz)/(15.*Delta);
	Dek[7][5] = 2.*(C0xy*C1xy + C0xy*C2xy + C0xz*C1xz + C0xz*C2xz + C0yz*C1yz + C0yz*C2yz + C1xy*C3xy + C1xz*C3xz + C1yz*C3yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(15.*Delta);
	Dek[7][6] = 2.*(C0xy*C0xy + C0xy*C2xy + C0xy*C3xy + C0xz*C0xz + C0xz*C2xz + C0xz*C3xz + C0yz*C0yz + C0yz*C2yz + C0yz*C3yz + 2.*C2xy*C3xy + 2.*C2xz*C3xz + 2.*C2yz*C3yz)/(15.*Delta);
	Dek[7][7] = 4.*(C0xy*C0xy + C0xy*C3xy + C0xz*C0xz + C0xz*C3xz + C0yz*C0yz + C0yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(15.*Delta);
	Dek[7][8] = 2.*(2*C0xy*C1xy + C0xy*C3xy + 2.*C0xz*C1xz + C0xz*C3xz + 2.*C0yz*C1yz + C0yz*C3yz + C1xy*C3xy + C1xz*C3xz + C1yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(15.*Delta);
	Dek[7][9] = 2.*(2*C0xy*C2xy + C0xy*C3xy + 2.*C0xz*C2xz + C0xz*C3xz + 2.*C0yz*C2yz + C0yz*C3yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(15.*Delta);
	Dek[8][0] = -(C0xy*C1xy + C0xy*C3xy + C0xz*C1xz + C0xz*C3xz + C0yz*C1yz + C0yz*C3yz)/(30.*Delta);
	Dek[8][1] = -(C1xy*C1xy - 3.*C1xy*C3xy + C1xz*C1xz - 3.*C1xz*C3xz + C1yz*C1yz - 3.*C1yz*C3yz)/(30.*Delta);
	Dek[8][2] = -(C1xy*C2xy + C1xz*C2xz + C1yz*C2yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(30.*Delta);
	Dek[8][3] = (3*C1xy*C3xy + 3.*C1xz*C3xz + 3.*C1yz*C3yz - C3xy*C3xy - C3xz*C3xz - C3yz*C3yz)/(30.*Delta);
	Dek[8][4] = 2.*(C0xy*C1xy + 2.*C0xy*C3xy + C0xz*C1xz + 2.*C0xz*C3xz + C0yz*C1yz + 2.*C0yz*C3yz + C1xy*C1xy + C1xy*C3xy + C1xz*C1xz + C1xz*C3xz + C1yz*C1yz + C1yz*C3yz)/(15.*Delta);
	Dek[8][5] = 2.*(C1xy*C1xy + C1xy*C2xy + C1xy*C3xy + C1xz*C1xz + C1xz*C2xz + C1xz*C3xz + C1yz*C1yz + C1yz*C2yz + C1yz*C3yz + 2.*C2xy*C3xy + 2.*C2xz*C3xz + 2.*C2yz*C3yz)/(15.*Delta);
	Dek[8][6] = 2.*(C0xy*C1xy + C0xy*C3xy + C0xz*C1xz + C0xz*C3xz + C0yz*C1yz + C0yz*C3yz + C1xy*C2xy + C1xz*C2xz + C1yz*C2yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz)/(15.*Delta);
	Dek[8][7] = 2.*(2*C0xy*C1xy + C0xy*C3xy + 2.*C0xz*C1xz + C0xz*C3xz + 2.*C0yz*C1yz + C0yz*C3yz + C1xy*C3xy + C1xz*C3xz + C1yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(15.*Delta);
	Dek[8][8] = 4.*(C1xy*C1xy + C1xy*C3xy + C1xz*C1xz + C1xz*C3xz + C1yz*C1yz + C1yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(15.*Delta);
	Dek[8][9] = 2.*(2*C1xy*C2xy + C1xy*C3xy + 2.*C1xz*C2xz + C1xz*C3xz + 2.*C1yz*C2yz + C1yz*C3yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(15.*Delta);
	Dek[9][0] = -(C0xy*C2xy + C0xy*C3xy + C0xz*C2xz + C0xz*C3xz + C0yz*C2yz + C0yz*C3yz)/(30.*Delta);
	Dek[9][1] = -(C1xy*C2xy + C1xy*C3xy + C1xz*C2xz + C1xz*C3xz + C1yz*C2yz + C1yz*C3yz)/(30.*Delta);
	Dek[9][2] = -(C2xy*C2xy - 3.*C2xy*C3xy + C2xz*C2xz - 3.*C2xz*C3xz + C2yz*C2yz - 3.*C2yz*C3yz)/(30.*Delta);
	Dek[9][3] =  (3*C2xy*C3xy + 3.*C2xz*C3xz + 3.*C2yz*C3yz - C3xy*C3xy - C3xz*C3xz - C3yz*C3yz)/(30.*Delta);
	Dek[9][4] = 2.*(C0xy*C2xy + C0xy*C3xy + C0xz*C2xz + C0xz*C3xz + C0yz*C2yz + C0yz*C3yz + C1xy*C2xy + C1xy*C3xy + C1xz*C2xz + C1xz*C3xz + C1yz*C2yz + C1yz*C3yz)/(15.*Delta);
	Dek[9][5] = 2.*(C1xy*C2xy + 2.*C1xy*C3xy + C1xz*C2xz + 2.*C1xz*C3xz + C1yz*C2yz + 2.*C1yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz)/(15.*Delta);
	Dek[9][6] = 2.*(C0xy*C2xy + 2.*C0xy*C3xy + C0xz*C2xz + 2.*C0xz*C3xz + C0yz*C2yz + 2.*C0yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz)/(15.*Delta);
	Dek[9][7] = 2.*(2*C0xy*C2xy + C0xy*C3xy + 2.*C0xz*C2xz + C0xz*C3xz + 2.*C0yz*C2yz + C0yz*C3yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(15.*Delta);
	Dek[9][8] = 2.*(2*C1xy*C2xy + C1xy*C3xy + 2.*C1xz*C2xz + C1xz*C3xz + 2.*C1yz*C2yz + C1yz*C3yz + C2xy*C3xy + C2xz*C3xz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(15.*Delta);
	Dek[9][9] = 4.*(C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(15.*Delta);

	return 0;
}

PetscErrorCode Calc3DP1b_Dek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **eP2_inter, PetscScalar **Dek)
{
	const PetscInt n0 = eP2_inter[en][0];
	const PetscInt n1 = eP2_inter[en][1];
	const PetscInt n2 = eP2_inter[en][2];
	const PetscInt n3 = eP2_inter[en][3];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
	const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];
	
	const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 + x1*y0*z2 
            -x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3  
            +x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0  
            +x3*y1*z2 + x3*y2*z0 - x3*y2*z1;

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

	Dek[0][0] = (945*C0xy*C0xy - 512*C0xy*C1xy + 945*C0xz*C0xz - 512*C0xz*C1xz + 945*C0yz*C0yz - 512*C0yz*C1yz + 512*C2xy*C2xy + 512*C2xy*C3xy + 512*C2xz*C2xz + 512*C2xz*C3xz + 512*C2yz*C2yz + 512*C2yz*C3yz + 512*C3xy*C3xy + 512*C3xz*C3xz + 512*C3yz*C3yz)/(5670*Delta);
	Dek[0][1] = (433*C0xy*C1xy + 433*C0xz*C1xz + 433*C0yz*C1yz + 512*C2xy*C2xy + 512*C2xy*C3xy + 512*C2xz*C2xz + 512*C2xz*C3xz + 512*C2yz*C2yz + 512*C2yz*C3yz + 512*C3xy*C3xy + 512*C3xz*C3xz + 512*C3yz*C3yz)/(5670*Delta);
	Dek[0][2] = (433*C0xy*C2xy + 433*C0xz*C2xz + 433*C0yz*C2yz + 512*C1xy*C1xy + 512*C1xy*C3xy + 512*C1xz*C1xz + 512*C1xz*C3xz + 512*C1yz*C1yz + 512*C1yz*C3yz + 512*C3xy*C3xy + 512*C3xz*C3xz + 512*C3yz*C3yz)/(5670*Delta);           
	Dek[0][3] = (433*C0xy*C3xy + 433*C0xz*C3xz + 433*C0yz*C3yz + 512*C1xy*C1xy + 512*C1xy*C2xy + 512*C1xz*C1xz + 512*C1xz*C2xz + 512*C1yz*C1yz + 512*C1yz*C2yz + 512*C2xy*C2xy + 512*C2xz*C2xz + 512*C2yz*C2yz)/(5670*Delta);           
	Dek[0][4] = -1024*(C0xy*C0xy + C0xy*C1xy + C0xy*C2xy + C0xy*C3xy + C0xz*C0xz + C0xz*C1xz + C0xz*C2xz + C0xz*C3xz + C0yz*C0yz + C0yz*C1yz + C0yz*C2yz + C0yz*C3yz + C1xy*C1xy + C1xy*C2xy + C1xy*C3xy + C1xz*C1xz + C1xz*C2xz + C1xz*C3xz + C1yz*C1yz + C1yz*C2yz + C1yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(2835*Delta);                                                                   

	Dek[1][0] =  (433*C0xy*C1xy + 433*C0xz*C1xz + 433*C0yz*C1yz + 512*C2xy*C2xy + 512*C2xy*C3xy + 512*C2xz*C2xz + 512*C2xz*C3xz + 512*C2yz*C2yz + 512*C2yz*C3yz + 512*C3xy*C3xy + 512*C3xz*C3xz + 512*C3yz*C3yz)/(5670*Delta);           
	Dek[1][1] = -(1008*C0xy*C1xy + 1008*C0xz*C1xz + 1008*C0yz*C1yz - 449*C1xy*C1xy + 496*C1xy*C2xy + 496*C1xy*C3xy - 449*C1xz*C1xz + 496*C1xz*C2xz + 496*C1xz*C3xz - 449*C1yz*C1yz + 496*C1yz*C2yz + 496*C1yz*C3yz - 512*C2xy*C2xy - 512*C2xy*C3xy - 512*C2xz*C2xz - 512*C2xz*C3xz - 512*C2yz*C2yz - 512*C2yz*C3yz - 512*C3xy*C3xy - 512*C3xz*C3xz - 512*C3yz*C3yz)/(5670*Delta);
	Dek[1][2] =  (512*C0xy*C0xy + 512*C0xy*C3xy + 512*C0xz*C0xz + 512*C0xz*C3xz + 512*C0yz*C0yz + 512*C0yz*C3yz + 433*C1xy*C2xy + 433*C1xz*C2xz + 433*C1yz*C2yz + 512*C3xy*C3xy + 512*C3xz*C3xz + 512*C3yz*C3yz)/(5670*Delta);           
	Dek[1][3] =  (512*C0xy*C0xy + 512*C0xy*C2xy + 512*C0xz*C0xz + 512*C0xz*C2xz + 512*C0yz*C0yz + 512*C0yz*C2yz + 433*C1xy*C3xy + 433*C1xz*C3xz + 433*C1yz*C3yz + 512*C2xy*C2xy + 512*C2xz*C2xz + 512*C2yz*C2yz)/(5670*Delta);          
	Dek[1][4] =  -1024*(C0xy*C0xy + C0xy*C1xy + C0xy*C2xy + C0xy*C3xy + C0xz*C0xz + C0xz*C1xz + C0xz*C2xz + C0xz*C3xz + C0yz*C0yz + C0yz*C1yz + C0yz*C2yz + C0yz*C3yz + C1xy*C1xy + C1xy*C2xy + C1xy*C3xy + C1xz*C1xz + C1xz*C2xz + C1xz*C3xz + C1yz*C1yz + C1yz*C2yz + C1yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(2835*Delta);                                                                   

	Dek[2][0] =  (433*C0xy*C2xy + 433*C0xz*C2xz + 433*C0yz*C2yz + 512*C1xy*C1xy + 512*C1xy*C3xy + 512*C1xz*C1xz + 512*C1xz*C3xz + 512*C1yz*C1yz + 512*C1yz*C3yz + 512*C3xy*C3xy + 512*C3xz*C3xz + 512*C3yz*C3yz)/(5670*Delta);           
	Dek[2][1] =  (512*C0xy*C0xy + 512*C0xy*C3xy + 512*C0xz*C0xz + 512*C0xz*C3xz + 512*C0yz*C0yz + 512*C0yz*C3yz + 433*C1xy*C2xy + 433*C1xz*C2xz + 433*C1yz*C2yz + 512*C3xy*C3xy + 512*C3xz*C3xz + 512*C3yz*C3yz)/(5670*Delta);           
	Dek[2][2] = -(1008*C0xy*C2xy + 1008*C0xz*C2xz + 1008*C0yz*C2yz - 512*C1xy*C1xy + 496*C1xy*C2xy - 512*C1xy*C3xy - 512*C1xz*C1xz + 496*C1xz*C2xz - 512*C1xz*C3xz - 512*C1yz*C1yz + 496*C1yz*C2yz - 512*C1yz*C3yz - 449*C2xy*C2xy + 496*C2xy*C3xy - 449*C2xz*C2xz + 496*C2xz*C3xz - 449*C2yz*C2yz + 496*C2yz*C3yz - 512*C3xy*C3xy - 512*C3xz*C3xz - 512*C3yz*C3yz)/(5670*Delta);
	Dek[2][3] =  (512*C0xy*C0xy + 512*C0xy*C1xy + 512*C0xz*C0xz + 512*C0xz*C1xz + 512*C0yz*C0yz + 512*C0yz*C1yz + 512*C1xy*C1xy + 512*C1xz*C1xz + 512*C1yz*C1yz + 433*C2xy*C3xy + 433*C2xz*C3xz + 433*C2yz*C3yz)/(5670*Delta);           
	Dek[2][4] =  -1024*(C0xy*C0xy + C0xy*C1xy + C0xy*C2xy + C0xy*C3xy + C0xz*C0xz + C0xz*C1xz + C0xz*C2xz + C0xz*C3xz + C0yz*C0yz + C0yz*C1yz + C0yz*C2yz + C0yz*C3yz + C1xy*C1xy + C1xy*C2xy + C1xy*C3xy + C1xz*C1xz + C1xz*C2xz + C1xz*C3xz + C1yz*C1yz + C1yz*C2yz + C1yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(2835*Delta);                                                                   

	Dek[3][0] =  (433*C0xy*C3xy + 433*C0xz*C3xz + 433*C0yz*C3yz + 512*C1xy*C1xy + 512*C1xy*C2xy + 512*C1xz*C1xz + 512*C1xz*C2xz + 512*C1yz*C1yz + 512*C1yz*C2yz + 512*C2xy*C2xy + 512*C2xz*C2xz + 512*C2yz*C2yz)/(5670*Delta);           
	Dek[3][1] =  (512*C0xy*C0xy + 512*C0xy*C2xy + 512*C0xz*C0xz + 512*C0xz*C2xz + 512*C0yz*C0yz + 512*C0yz*C2yz + 433*C1xy*C3xy + 433*C1xz*C3xz + 433*C1yz*C3yz + 512*C2xy*C2xy + 512*C2xz*C2xz + 512*C2yz*C2yz)/(5670*Delta);           
	Dek[3][2] =  (512*C0xy*C0xy + 512*C0xy*C1xy + 512*C0xz*C0xz + 512*C0xz*C1xz + 512*C0yz*C0yz + 512*C0yz*C1yz + 512*C1xy*C1xy + 512*C1xz*C1xz + 512*C1yz*C1yz + 433*C2xy*C3xy + 433*C2xz*C3xz + 433*C2yz*C3yz)/(5670*Delta);           
	Dek[3][3] =  -(1008*C0xy*C3xy + 1008*C0xz*C3xz + 1008*C0yz*C3yz - 512*C1xy*C1xy - 512*C1xy*C2xy + 496*C1xy*C3xy - 512*C1xz*C1xz - 512*C1xz*C2xz + 496*C1xz*C3xz - 512*C1yz*C1yz - 512*C1yz*C2yz + 496*C1yz*C3yz - 512*C2xy*C2xy + 496*C2xy*C3xy - 512*C2xz*C2xz + 496*C2xz*C3xz - 512*C2yz*C2yz + 496*C2yz*C3yz - 449*C3xy*C3xy - 449*C3xz*C3xz - 449*C3yz*C3yz)/(5670*Delta);
	Dek[3][4] =  -1024*(C0xy*C0xy + C0xy*C1xy + C0xy*C2xy + C0xy*C3xy + C0xz*C0xz + C0xz*C1xz + C0xz*C2xz + C0xz*C3xz + C0yz*C0yz + C0yz*C1yz + C0yz*C2yz + C0yz*C3yz + C1xy*C1xy + C1xy*C2xy + C1xy*C3xy + C1xz*C1xz + C1xz*C2xz + C1xz*C3xz + C1yz*C1yz + C1yz*C2yz + C1yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(2835*Delta);

	Dek[4][0] =  -1024*(C0xy*C0xy + C0xy*C1xy + C0xy*C2xy + C0xy*C3xy + C0xz*C0xz + C0xz*C1xz + C0xz*C2xz + C0xz*C3xz + C0yz*C0yz + C0yz*C1yz + C0yz*C2yz + C0yz*C3yz + C1xy*C1xy + C1xy*C2xy + C1xy*C3xy + C1xz*C1xz + C1xz*C2xz + C1xz*C3xz + C1yz*C1yz + C1yz*C2yz + C1yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(2835*Delta);                                                                   
	Dek[4][1] =  -1024*(C0xy*C0xy + C0xy*C1xy + C0xy*C2xy + C0xy*C3xy + C0xz*C0xz + C0xz*C1xz + C0xz*C2xz + C0xz*C3xz + C0yz*C0yz + C0yz*C1yz + C0yz*C2yz + C0yz*C3yz + C1xy*C1xy + C1xy*C2xy + C1xy*C3xy + C1xz*C1xz + C1xz*C2xz + C1xz*C3xz + C1yz*C1yz + C1yz*C2yz + C1yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(2835*Delta);                                                                   
	Dek[4][2] =  -1024*(C0xy*C0xy + C0xy*C1xy + C0xy*C2xy + C0xy*C3xy + C0xz*C0xz + C0xz*C1xz + C0xz*C2xz + C0xz*C3xz + C0yz*C0yz + C0yz*C1yz + C0yz*C2yz + C0yz*C3yz + C1xy*C1xy + C1xy*C2xy + C1xy*C3xy + C1xz*C1xz + C1xz*C2xz + C1xz*C3xz + C1yz*C1yz + C1yz*C2yz + C1yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(2835*Delta);                                                                   
	Dek[4][3] =  -1024*(C0xy*C0xy + C0xy*C1xy + C0xy*C2xy + C0xy*C3xy + C0xz*C0xz + C0xz*C1xz + C0xz*C2xz + C0xz*C3xz + C0yz*C0yz + C0yz*C1yz + C0yz*C2yz + C0yz*C3yz + C1xy*C1xy + C1xy*C2xy + C1xy*C3xy + C1xz*C1xz + C1xz*C2xz + C1xz*C3xz + C1yz*C1yz + C1yz*C2yz + C1yz*C3yz + C2xy*C2xy + C2xy*C3xy + C2xz*C2xz + C2xz*C3xz + C2yz*C2yz + C2yz*C3yz + C3xy*C3xy + C3xz*C3xz + C3yz*C3yz)/(2835*Delta);                                                                   
	Dek[4][4] =  -4096*(C0xy*C1xy + C0xz*C1xz + C0yz*C1yz - C2xy*C2xy - C2xy*C3xy - C2xz*C2xz - C2xz*C3xz - C2yz*C2yz - C2yz*C3yz - C3xy*C3xy - C3xz*C3xz - C3yz*C3yz)/(2835*Delta);


	return 0;
}
#endif