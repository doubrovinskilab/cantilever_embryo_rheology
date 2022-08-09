#ifndef ElementalVectorSource_C
#define ElementalVectorSource_C

PetscErrorCode Calc2DP1_Fek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
 PetscInt **e_inter, PetscScalar *S, PetscScalar *Fek, const PetscBool TypeS)
{
	const PetscInt n0 = e_inter[en][0];
	const PetscInt n1 = e_inter[en][1];
	const PetscInt n2 = e_inter[en][2];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2];
	
	const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;

	PetscScalar f0, f1, f2;
	if (TypeS) {
		f0 = S[n0];
		f1 = S[n1];
		f2 = S[n2]; 	}
	else 	{
		f0 = S[0];
		f1 = S[1];
		f2 = S[2];	}

	Fek[0] = Delta*(2.*f0 + f1 + f2)/24.;
  Fek[1] = Delta*(f0 + 2.*f1 + f2)/24.;
  Fek[2] = Delta*(f0 + f1 + 2.*f2)/24.;

  return 0;
}

PetscErrorCode Calc2DP2_Fek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
 PetscInt **e_inter, PetscScalar *S, PetscScalar *Fek, const PetscBool TypeS)
{
	const PetscInt n0 = e_inter[en][0];
	const PetscInt n1 = e_inter[en][1];
	const PetscInt n2 = e_inter[en][2];
	const PetscInt n3 = e_inter[en][3];
	const PetscInt n4 = e_inter[en][4];
	const PetscInt n5 = e_inter[en][5];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2];

	const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;

	PetscScalar f0, f1, f2, f01, f12, f02;
	if (TypeS) {
		f0  = S[n0];
		f1  = S[n1];
		f2  = S[n2];
		f01 = S[n3];
		f12 = S[n4];
		f02 = S[n5]; }
	else {
		f0  = S[0];
		f1  = S[1];
		f2  = S[2];
		f01 = S[3];
		f12 = S[4];
		f02 = S[5]; }

	Fek[0] =  Delta*(6.*f0 - f1 - 4.*f12 - f2)/360.;
	Fek[1] = -Delta*(f0 + 4.*f02 - 6.*f1 + f2)/360.;
	Fek[2] = -Delta*(f0 + 4.*f01 + f1 - 6.*f2)/360.;
	Fek[3] =  Delta*(8.*f01 + 4.*f02 + 4.*f12 - f2)/90.;
	Fek[4] = -Delta*(f0 - 4.*f01 - 4.*f02 - 8.*f12)/90.;
	Fek[5] =  Delta*(4.*f01 + 8.*f02 - f1 + 4.*f12)/90.;

	return 0;
}

PetscErrorCode Calc2DP1b_Fek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
 PetscInt **e_inter, PetscScalar *S, PetscScalar *Fek, const PetscBool TypeS)
{

	const PetscInt n0   = e_inter[en][0];
	const PetscInt n1   = e_inter[en][1];
	const PetscInt n2   = e_inter[en][2];
	const PetscInt n012 = e_inter[en][3];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2];
	
	const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;

	PetscScalar f0, f1, f2, f012;
	if (TypeS) {
		f0 = S[n0];
		f1 = S[n1];
		f2 = S[n2];
		f012 = S[n012]; }
	else {
		f0 = S[0];
		f1 = S[1];
		f2 = S[2];
		f012 = S[3]; }

  Fek[0] =    Delta*(83.*f0 + 45.*f012 + 13.*f1 + 13.*f2)/1680.;
  Fek[1] =    Delta*(13.*f0 + 45.*f012 + 83.*f1 + 13.*f2)/1680.;
  Fek[2] =    Delta*(13.*f0 + 45.*f012 + 13.*f1 + 83.*f2)/1680.;
  Fek[3] = 3.*Delta*(5.*f0  + 27.*f012 +  5.*f1 +  5.*f2)/560.;

  return 0;
}


PetscErrorCode Calc3DP1_Fek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
 PetscInt **e_inter, PetscScalar *S, PetscScalar *Fek, const PetscBool TypeS)
{

	const PetscInt n0 = e_inter[en][0];
	const PetscInt n1 = e_inter[en][1];
	const PetscInt n2 = e_inter[en][2];
	const PetscInt n3 = e_inter[en][3];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
	const PetscScalar x3 = n_x[n2], y3 = n_y[n2], z3 = n_z[n3];

	const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 + x1*y0*z2 - x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3 + x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0 + x3*y1*z2 + x3*y2*z0 - x3*y2*z1;

	PetscScalar f0, f1, f2, f3;
	if (TypeS) {
		f0 = S[n0];
		f1 = S[n1];
		f2 = S[n2];
		f3 = S[n3];}
	else {
		f0 = S[0];
		f1 = S[1];
		f2 = S[2];
		f3 = S[3];}

	Fek[0] = Delta*(2.*f0 + f1 + f2 + f3)/120.;
  Fek[1] = Delta*(f0 + 2.*f1 + f2 + f3)/120.;
  Fek[2] = Delta*(f0 + f1 + 2.*f2 + f3)/120.;
  Fek[3] = Delta*(f0 + f1 + f2 + 2.*f3)/120.;

	return 0;
}

PetscErrorCode Calc3DP2_Fek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
 PetscInt **e_inter, PetscScalar *S, PetscScalar *Fek, const PetscBool TypeS)
{
	const PetscInt n0 = e_inter[en][0];
	const PetscInt n1 = e_inter[en][1];
	const PetscInt n2 = e_inter[en][2];
	const PetscInt n3 = e_inter[en][3];
	const PetscInt n4 = e_inter[en][4];
	const PetscInt n5 = e_inter[en][5];
	const PetscInt n6 = e_inter[en][6];
	const PetscInt n7 = e_inter[en][7];
	const PetscInt n8 = e_inter[en][8];
	const PetscInt n9 = e_inter[en][9];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
	const PetscScalar x3 = n_x[n2], y3 = n_y[n2], z3 = n_z[n3];

	const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 
            + x1*y0*z2 -x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 
            + x2*y0*z3 +x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 
            - x3*y1*z0 +x3*y1*z2 + x3*y2*z0 - x3*y2*z1;

	PetscScalar f0, f1, f2, f3, f01, f12, f02, f03, f13, f23;
	if (TypeS) {
		f0  = S[n0];
		f1  = S[n1];
		f2  = S[n2];
		f3  = S[n3];
		f01 = S[n4];
		f12 = S[n5];
		f02 = S[n6];
		f03 = S[n7];
		f13 = S[n8];
		f23 = S[n9];}
	else {
		f0  = S[0];
		f1  = S[1];
		f2  = S[2];
		f3  = S[3];
		f01 = S[4];
		f12 = S[5];
		f02 = S[6];
		f03 = S[7];
		f13 = S[8];
		f23 = S[9];}

	Fek[0] =  Delta*(6*f0 - 4.*f01 - 4.*f02 - 4.*f03 + f1 - 6.*f12 - 6.*f13 + f2 - 6.*f23 + f3)/2520.;
	Fek[1] =  Delta*(f0 - 4.*f01 - 6.*f02 - 6.*f03 + 6.*f1 - 4.*f12 - 4.*f13 + f2 - 6.*f23 + f3)/2520.;
	Fek[2] =  Delta*(f0 - 6.*f01 - 4.*f02 - 6.*f03 + f1 - 4.*f12 - 6.*f13 + 6.*f2 - 4.*f23 + f3)/2520.;
	Fek[3] =  Delta*(f0 - 6.*f01 - 6.*f02 - 4.*f03 + f1 - 6.*f12 - 4.*f13 + f2 - 4.*f23 + 6.*f3)/2520.;
	Fek[4] = -Delta*(2.*f0 - 16.*f01 - 8.*f02 - 8.*f03 + 2.*f1 - 8.*f12 - 8.*f13 + 3.*f2 - 4.*f23 + 3.*f3)/1260.;
	Fek[5] = -Delta*(3.*f0 - 8.*f01 - 8.*f02 - 4.*f03 + 2.*f1 - 16.*f12 - 8.*f13 + 2.*f2 - 8.*f23 + 3.*f3)/1260.;
	Fek[6] = -Delta*(2.*f0 - 8.*f01 - 16.*f02 - 8.*f03 + 3.*f1 - 8.*f12 - 4.*f13 + 2.*f2 - 8.*f23 + 3.*f3)/1260.;
	Fek[7] = -Delta*(2.*f0 - 8.*f01 - 8.*f02 - 16.*f03 + 3.*f1 - 4.*f12 - 8.*f13 + 3.*f2 - 8.*f23 + 2.*f3)/1260.;
	Fek[8] = -Delta*(3.*f0 - 8.*f01 - 4.*f02 - 8.*f03 + 2.*f1 - 8.*f12 - 16.*f13 + 3.*f2 - 8.*f23 + 2.*f3)/1260.;
	Fek[9] = -Delta*(3.*f0 - 4.*f01 - 8.*f02 - 8.*f03 + 3.*f1 - 8.*f12 - 8.*f13 + 2.*f2 - 16.*f23 + 2.*f3)/1260.;

	return 0;
}


PetscErrorCode Calc3DP1b_Fek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
 PetscInt **e_inter, PetscScalar *S, PetscScalar *Fek, const PetscBool TypeS)
{
	const PetscInt n0 = e_inter[en][0];
	const PetscInt n1 = e_inter[en][1];
	const PetscInt n2 = e_inter[en][2];
	const PetscInt n3 = e_inter[en][3];
	const PetscInt n4 = e_inter[en][4];

	const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
	const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
	const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
	const PetscScalar x3 = n_x[n2], y3 = n_y[n2], z3 = n_z[n3];

	const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 
            + x1*y0*z2 -x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 
            + x2*y0*z3 +x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 
            - x3*y1*z0 +x3*y1*z2 + x3*y2*z0 - x3*y2*z1;

	PetscScalar f0, f1, f2, f3, f0123;
	if (TypeS) {
		f0  	= S[n0];
		f1  	= S[n1];
		f2  	= S[n2];
		f3  	= S[n3];
		f0123 = S[n4];}
	else {
		f0  	= S[0];
		f1  	= S[1];
		f2  	= S[2];
		f3  	= S[3];
		f0123 = S[4];}

  Fek[0] =  Delta*(14918*f0 + 7648*f0123 + 4523*f1 + 4523*f2 + 4523*f3)/1247400.;
	Fek[1] =  Delta*(4523*f0 + 7648*f0123 + 14918*f1 + 4523*f2 + 4523*f3)/1247400.;
	Fek[2] =  Delta*(4523*f0 + 7648*f0123 + 4523*f1 + 14918*f2 + 4523*f3)/1247400.;
	Fek[3] =  Delta*(4523*f0 + 7648*f0123 + 4523*f1 + 4523*f2 + 14918*f3)/1247400.;
	Fek[4] =  4*Delta*(239*f0 + 1024*f0123 + 239*f1 + 239*f2 + 239*f3)/155925.;

	return 0;
}
#endif