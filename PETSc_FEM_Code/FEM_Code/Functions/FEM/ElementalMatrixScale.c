#ifndef ElementalMatrixScale_C
#define ElementalMatrixScale_C

PetscErrorCode Calc2DP1P1_Scale(PetscScalar **Sek, const PetscScalar S)
{
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;
  return 0;
}

PetscErrorCode Calc2DP2P2_Scale(PetscScalar **Sek, const PetscScalar S)
{
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;  Sek[0][4] = Sek[0][4]*S;  Sek[0][5] = Sek[0][5]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;  Sek[1][4] = Sek[1][4]*S;  Sek[1][5] = Sek[1][5]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;  Sek[2][4] = Sek[2][4]*S;  Sek[2][5] = Sek[2][5]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;  Sek[3][3] = Sek[3][3]*S;  Sek[3][4] = Sek[3][4]*S;  Sek[3][5] = Sek[3][5]*S;
  Sek[4][0] = Sek[4][0]*S;  Sek[4][1] = Sek[4][1]*S;  Sek[4][2] = Sek[4][2]*S;  Sek[4][3] = Sek[4][3]*S;  Sek[4][4] = Sek[4][4]*S;  Sek[4][5] = Sek[4][5]*S;  
  Sek[5][0] = Sek[5][0]*S;  Sek[5][1] = Sek[5][1]*S;  Sek[5][2] = Sek[5][2]*S;  Sek[5][3] = Sek[5][3]*S;  Sek[5][4] = Sek[5][4]*S;  Sek[5][5] = Sek[5][5]*S;
	return 0;
}

PetscErrorCode Calc2DP2P1_Scale(PetscScalar **Sek, const PetscScalar S)
{
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;
  Sek[4][0] = Sek[4][0]*S;  Sek[4][1] = Sek[4][1]*S;  Sek[4][2] = Sek[4][2]*S;
  Sek[5][0] = Sek[5][0]*S;  Sek[5][1] = Sek[5][1]*S;  Sek[5][2] = Sek[5][2]*S;
	return 0;
}
PetscErrorCode Calc2DP1P2_Scale(PetscScalar **Sek, const PetscScalar S)
{
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;  Sek[0][4] = Sek[0][4]*S;  Sek[0][5] = Sek[0][5]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;  Sek[1][4] = Sek[1][4]*S;  Sek[1][5] = Sek[1][5]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;  Sek[2][4] = Sek[2][4]*S;  Sek[2][5] = Sek[2][5]*S;
	return 0;
}


PetscErrorCode Calc2DP1bP1_Scale(PetscScalar **Sek, const PetscScalar S)
{
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;
  return 0;
}

PetscErrorCode Calc2DP1P1b_Scale(PetscScalar **Sek, const PetscScalar S)
{ 
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;
  return 0;
}

PetscErrorCode Calc2DP1bP1b_Scale(PetscScalar **Sek, const PetscScalar S)
{ 
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;  Sek[3][3] = Sek[3][3]*S;
  return 0;
}

PetscErrorCode Calc3DP1P1_Scale(PetscScalar **Sek, const PetscScalar S)
{ 
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;  Sek[3][3] = Sek[3][3]*S;
  return 0;
}PetscErrorCode Calc3DP2P2_Scale(PetscScalar **Sek, const PetscScalar S)
{
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;  Sek[0][4] = Sek[0][4]*S;
  Sek[0][5] = Sek[0][5]*S;  Sek[0][6] = Sek[0][6]*S;  Sek[0][7] = Sek[0][7]*S;  Sek[0][8] = Sek[0][8]*S;  Sek[0][9] = Sek[0][9]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;  Sek[1][4] = Sek[1][4]*S;
  Sek[1][5] = Sek[1][5]*S;  Sek[1][6] = Sek[1][6]*S;  Sek[1][7] = Sek[1][7]*S;  Sek[1][8] = Sek[1][8]*S;  Sek[1][9] = Sek[1][9]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;  Sek[2][4] = Sek[2][4]*S;
  Sek[2][5] = Sek[2][5]*S;  Sek[2][6] = Sek[2][6]*S;  Sek[2][7] = Sek[2][7]*S;  Sek[2][8] = Sek[2][8]*S;  Sek[2][9] = Sek[2][9]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;  Sek[3][3] = Sek[3][3]*S;  Sek[3][4] = Sek[3][4]*S;
  Sek[3][5] = Sek[3][5]*S;  Sek[3][6] = Sek[3][6]*S;  Sek[3][7] = Sek[3][7]*S;  Sek[3][8] = Sek[3][8]*S;  Sek[3][9] = Sek[3][9]*S;
  Sek[4][0] = Sek[4][0]*S;  Sek[4][1] = Sek[4][1]*S;  Sek[4][2] = Sek[4][2]*S;  Sek[4][3] = Sek[4][3]*S;  Sek[4][4] = Sek[4][4]*S;
  Sek[4][5] = Sek[4][5]*S;  Sek[4][6] = Sek[4][6]*S;  Sek[4][7] = Sek[4][7]*S;  Sek[4][8] = Sek[4][8]*S;  Sek[4][9] = Sek[4][9]*S;
  Sek[5][0] = Sek[5][0]*S;  Sek[5][1] = Sek[5][1]*S;  Sek[5][2] = Sek[5][2]*S;  Sek[5][3] = Sek[5][3]*S;  Sek[5][4] = Sek[5][4]*S;
  Sek[5][5] = Sek[5][5]*S;  Sek[5][6] = Sek[5][6]*S;  Sek[5][7] = Sek[5][7]*S;  Sek[5][8] = Sek[5][8]*S;  Sek[5][9] = Sek[5][9]*S;
  Sek[6][0] = Sek[6][0]*S;  Sek[6][1] = Sek[6][1]*S;  Sek[6][2] = Sek[6][2]*S;  Sek[6][3] = Sek[6][3]*S;  Sek[6][4] = Sek[6][4]*S;
  Sek[6][5] = Sek[6][5]*S;  Sek[6][6] = Sek[6][6]*S;  Sek[6][7] = Sek[6][7]*S;  Sek[6][8] = Sek[6][8]*S;  Sek[6][9] = Sek[6][9]*S;
  Sek[7][0] = Sek[7][0]*S;  Sek[7][1] = Sek[7][1]*S;  Sek[7][2] = Sek[7][2]*S;  Sek[7][3] = Sek[7][3]*S;  Sek[7][4] = Sek[7][4]*S;
  Sek[7][5] = Sek[7][5]*S;  Sek[7][6] = Sek[7][6]*S;  Sek[7][7] = Sek[7][7]*S;  Sek[7][8] = Sek[7][8]*S;  Sek[7][9] = Sek[7][9]*S;
  Sek[8][0] = Sek[8][0]*S;  Sek[8][1] = Sek[8][1]*S;  Sek[8][2] = Sek[8][2]*S;  Sek[8][3] = Sek[8][3]*S;  Sek[8][4] = Sek[8][4]*S;
  Sek[8][5] = Sek[8][5]*S;  Sek[8][6] = Sek[8][6]*S;  Sek[8][7] = Sek[8][7]*S;  Sek[8][8] = Sek[8][8]*S;  Sek[8][9] = Sek[8][9]*S;
  Sek[9][0] = Sek[9][0]*S;  Sek[9][1] = Sek[9][1]*S;  Sek[9][2] = Sek[9][2]*S;  Sek[9][3] = Sek[9][3]*S;  Sek[9][4] = Sek[9][4]*S;
  Sek[9][5] = Sek[9][5]*S;  Sek[9][6] = Sek[9][6]*S;  Sek[9][7] = Sek[9][7]*S;  Sek[9][8] = Sek[9][8]*S;  Sek[9][9] = Sek[9][9]*S;

	return 0;
}
PetscErrorCode Calc3DP2P1_Scale(PetscScalar **Sek, const PetscScalar S)
{

  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;  Sek[3][3] = Sek[3][3]*S;
  Sek[4][0] = Sek[4][0]*S;  Sek[4][1] = Sek[4][1]*S;  Sek[4][2] = Sek[4][2]*S;  Sek[4][3] = Sek[4][3]*S;
  Sek[5][0] = Sek[5][0]*S;  Sek[5][1] = Sek[5][1]*S;  Sek[5][2] = Sek[5][2]*S;  Sek[5][3] = Sek[5][3]*S;
  Sek[6][0] = Sek[6][0]*S;  Sek[6][1] = Sek[6][1]*S;  Sek[6][2] = Sek[6][2]*S;  Sek[6][3] = Sek[6][3]*S;
  Sek[7][0] = Sek[7][0]*S;  Sek[7][1] = Sek[7][1]*S;  Sek[7][2] = Sek[7][2]*S;  Sek[7][3] = Sek[7][3]*S;
  Sek[8][0] = Sek[8][0]*S;  Sek[8][1] = Sek[8][1]*S;  Sek[8][2] = Sek[8][2]*S;  Sek[8][3] = Sek[8][3]*S;
  Sek[9][0] = Sek[9][0]*S;  Sek[9][1] = Sek[9][1]*S;  Sek[9][2] = Sek[9][2]*S;  Sek[9][3] = Sek[9][3]*S;
 	return 0;
}

PetscErrorCode Calc3DP1P2_Scale(PetscScalar **Sek, const PetscScalar S)
{
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;  Sek[0][4] = Sek[0][4]*S;
  Sek[0][5] = Sek[0][5]*S;  Sek[0][6] = Sek[0][6]*S;  Sek[0][7] = Sek[0][7]*S;  Sek[0][8] = Sek[0][8]*S;  Sek[0][9] = Sek[0][9]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;  Sek[1][4] = Sek[1][4]*S;
  Sek[1][5] = Sek[1][5]*S;  Sek[1][6] = Sek[1][6]*S;  Sek[1][7] = Sek[1][7]*S;  Sek[1][8] = Sek[1][8]*S;  Sek[1][9] = Sek[1][9]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;  Sek[2][4] = Sek[2][4]*S;
  Sek[2][5] = Sek[2][5]*S;  Sek[2][6] = Sek[2][6]*S;  Sek[2][7] = Sek[2][7]*S;  Sek[2][8] = Sek[2][8]*S;  Sek[2][9] = Sek[2][9]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;  Sek[3][3] = Sek[3][3]*S;  Sek[3][4] = Sek[3][4]*S;
  Sek[3][5] = Sek[3][5]*S;  Sek[3][6] = Sek[3][6]*S;  Sek[3][7] = Sek[3][7]*S;  Sek[3][8] = Sek[3][8]*S;  Sek[3][9] = Sek[3][9]*S;
	return 0;
}PetscErrorCode Calc3DP1bP1_Scale(PetscScalar **Sek, const PetscScalar S)
{
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;  Sek[3][3] = Sek[3][3]*S;
  Sek[4][0] = Sek[4][0]*S;  Sek[4][1] = Sek[4][1]*S;  Sek[4][2] = Sek[4][2]*S;  Sek[4][3] = Sek[4][3]*S;
  return 0;
}
PetscErrorCode Calc3DP1P1b_Scale(PetscScalar **Sek, const PetscScalar S)
{
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;  Sek[0][4] = Sek[0][4]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;  Sek[1][4] = Sek[1][4]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;  Sek[2][4] = Sek[2][4]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;  Sek[3][3] = Sek[3][3]*S;  Sek[3][4] = Sek[3][4]*S;
  return 0;
}

PetscErrorCode Calc3DP1bP1b_Scale(PetscScalar **Sek, const PetscScalar S)
{
  Sek[0][0] = Sek[0][0]*S;  Sek[0][1] = Sek[0][1]*S;  Sek[0][2] = Sek[0][2]*S;  Sek[0][3] = Sek[0][3]*S;  Sek[0][4] = Sek[0][4]*S;
  Sek[1][0] = Sek[1][0]*S;  Sek[1][1] = Sek[1][1]*S;  Sek[1][2] = Sek[1][2]*S;  Sek[1][3] = Sek[1][3]*S;  Sek[1][4] = Sek[1][4]*S;
  Sek[2][0] = Sek[2][0]*S;  Sek[2][1] = Sek[2][1]*S;  Sek[2][2] = Sek[2][2]*S;  Sek[2][3] = Sek[2][3]*S;  Sek[2][4] = Sek[2][4]*S;
  Sek[3][0] = Sek[3][0]*S;  Sek[3][1] = Sek[3][1]*S;  Sek[3][2] = Sek[3][2]*S;  Sek[3][3] = Sek[3][3]*S;  Sek[3][4] = Sek[3][4]*S;
  Sek[4][0] = Sek[4][0]*S;  Sek[4][1] = Sek[4][1]*S;  Sek[4][2] = Sek[4][2]*S;  Sek[4][3] = Sek[4][3]*S;  Sek[4][4] = Sek[4][4]*S;
  return 0;
}

#endif