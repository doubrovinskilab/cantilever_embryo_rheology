#ifndef AddVariablesToGlobalStokesFullLinearSystem_H
#define AddVariablesToGlobalStokesFullLinearSystem_H

PetscErrorCode AddMomEqs2D(Mat *A, Vec *b, const PetscInt U_S_inter, const PetscInt P_S_inter, const PetscInt n, 
  PetscInt **GI_V_elem, PetscScalar **Dek, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z, PetscScalar *Fek_x, 
  PetscScalar *Fek_y, PetscScalar *Fek_z)
{
  MatSetValues(*A, 1 , &GI_V_elem[0][n], U_S_inter, GI_V_elem[0], Dek[n],   ADD_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[0][n], P_S_inter, GI_V_elem[2], Gek_x[n], ADD_VALUES);  
  VecSetValue(*b, GI_V_elem[0][n], Fek_x[n], ADD_VALUES);                                 

  MatSetValues(*A, 1 , &GI_V_elem[1][n], U_S_inter, GI_V_elem[1], Dek[n],   ADD_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[1][n], P_S_inter, GI_V_elem[2], Gek_y[n], ADD_VALUES);  
  VecSetValue(*b, GI_V_elem[1][n], Fek_y[n], ADD_VALUES);                                 
  return 0;
}

PetscErrorCode AddMomEqs3D(Mat *A, Vec *b, const PetscInt U_S_inter, const PetscInt P_S_inter, const PetscInt n, 
  PetscInt **GI_V_elem, PetscScalar **Dek, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z, PetscScalar *Fek_x, 
  PetscScalar *Fek_y, PetscScalar *Fek_z)
{
  MatSetValues(*A, 1 , &GI_V_elem[0][n], U_S_inter, GI_V_elem[0], Dek[n],   ADD_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[0][n], P_S_inter, GI_V_elem[3], Gek_x[n], ADD_VALUES);  
  VecSetValue(*b, GI_V_elem[0][n], Fek_x[n], ADD_VALUES);                                 

  MatSetValues(*A, 1 , &GI_V_elem[1][n], U_S_inter, GI_V_elem[1], Dek[n],   ADD_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[1][n], P_S_inter, GI_V_elem[3], Gek_y[n], ADD_VALUES);  
  VecSetValue(*b, GI_V_elem[1][n], Fek_y[n], ADD_VALUES);                                 

  MatSetValues(*A, 1 , &GI_V_elem[2][n], U_S_inter, GI_V_elem[2], Dek[n],   ADD_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[2][n], P_S_inter, GI_V_elem[3], Gek_z[n], ADD_VALUES);  
  VecSetValue(*b, GI_V_elem[2][n], Fek_z[n], ADD_VALUES);                                 
  return 0;
}

PetscErrorCode AddContEqs2DNoPress(Mat *A, const PetscInt U_S_inter, const PetscInt P_S_inter, const PetscInt n, 
  PetscInt **GI_V_elem, PetscScalar **Mpek, PetscScalar **Gek_xT, PetscScalar **Gek_yT, PetscScalar **Gek_zT)
{
  MatSetValues(*A, 1 , &GI_V_elem[2][n], U_S_inter, GI_V_elem[0], Gek_xT[n], ADD_VALUES);
  MatSetValues(*A, 1 , &GI_V_elem[2][n], U_S_inter, GI_V_elem[1], Gek_yT[n], ADD_VALUES);
  return 0;
}

PetscErrorCode AddContEqs3DNoPress(Mat *A, const PetscInt U_S_inter, const PetscInt P_S_inter, const PetscInt n, 
  PetscInt **GI_V_elem, PetscScalar **Mpek, PetscScalar **Gek_xT, PetscScalar **Gek_yT, PetscScalar **Gek_zT)
{  
  MatSetValues(*A, 1 , &GI_V_elem[3][n], U_S_inter, GI_V_elem[0], Gek_xT[n], ADD_VALUES);
  MatSetValues(*A, 1 , &GI_V_elem[3][n], U_S_inter, GI_V_elem[1], Gek_yT[n], ADD_VALUES);
  MatSetValues(*A, 1 , &GI_V_elem[3][n], U_S_inter, GI_V_elem[2], Gek_zT[n], ADD_VALUES);
  return 0;
}

PetscErrorCode AddContEqs2DZeroDiagPress(Mat *A, const PetscInt U_S_inter, const PetscInt P_S_inter, const PetscInt n, 
  PetscInt **GI_V_elem, PetscScalar **Mpek, PetscScalar **Gek_xT, PetscScalar **Gek_yT, PetscScalar **Gek_zT)
{
  PetscScalar zero = 0.0;
  MatSetValues(*A, 1 , &GI_V_elem[2][n], U_S_inter, GI_V_elem[0], Gek_xT[n], ADD_VALUES); 
  MatSetValues(*A, 1 , &GI_V_elem[2][n], U_S_inter, GI_V_elem[1], Gek_yT[n], ADD_VALUES); 
  MatSetValues(*A, 1 , &GI_V_elem[2][n],         1,&GI_V_elem[2][n],  &zero, ADD_VALUES); 
  return 0;
}

PetscErrorCode AddContEqs3DZeroDiagPress(Mat *A, const PetscInt U_S_inter, const PetscInt P_S_inter, const PetscInt n, 
  PetscInt **GI_V_elem, PetscScalar **Mpek, PetscScalar **Gek_xT, PetscScalar **Gek_yT, PetscScalar **Gek_zT)
{
  PetscScalar zero = 0.0;
  MatSetValues(*A, 1 , &GI_V_elem[3][n], U_S_inter, GI_V_elem[0],   Gek_xT[n], ADD_VALUES); 
  MatSetValues(*A, 1 , &GI_V_elem[3][n], U_S_inter, GI_V_elem[1],   Gek_yT[n], ADD_VALUES); 
  MatSetValues(*A, 1 , &GI_V_elem[3][n], U_S_inter, GI_V_elem[2],   Gek_zT[n], ADD_VALUES); 
  MatSetValues(*A, 1 , &GI_V_elem[3][n],         1,&GI_V_elem[3][n],    &zero, ADD_VALUES); 
  return 0;
}

PetscErrorCode AddContEqs2DMassPress(Mat *A, const PetscInt U_S_inter, const PetscInt P_S_inter, const PetscInt n, 
  PetscInt **GI_V_elem, PetscScalar **Mpek, PetscScalar **Gek_xT, PetscScalar **Gek_yT, PetscScalar **Gek_zT)
{
  MatSetValues(*A, 1 , &GI_V_elem[2][n], U_S_inter, GI_V_elem[0], Gek_xT[n], ADD_VALUES); 
  MatSetValues(*A, 1 , &GI_V_elem[2][n], U_S_inter, GI_V_elem[1], Gek_yT[n], ADD_VALUES); 
  MatSetValues(*A, 1 , &GI_V_elem[2][n], P_S_inter, GI_V_elem[2], Mpek[n]  , ADD_VALUES); 
  return 0;
}

PetscErrorCode AddContEqs3DMassPress(Mat *A, const PetscInt U_S_inter, const PetscInt P_S_inter, const PetscInt n, 
  PetscInt **GI_V_elem, PetscScalar **Mpek, PetscScalar **Gek_xT, PetscScalar **Gek_yT, PetscScalar **Gek_zT)
{
  MatSetValues(*A, 1 , &GI_V_elem[3][n], U_S_inter, GI_V_elem[0],   Gek_xT[n], ADD_VALUES); 
  MatSetValues(*A, 1 , &GI_V_elem[3][n], U_S_inter, GI_V_elem[1],   Gek_yT[n], ADD_VALUES); 
  MatSetValues(*A, 1 , &GI_V_elem[3][n], U_S_inter, GI_V_elem[2],   Gek_zT[n], ADD_VALUES); 
  MatSetValues(*A, 1 , &GI_V_elem[3][n], P_S_inter, GI_V_elem[3],   Mpek[n]  , ADD_VALUES); 
  return 0;
}

PetscErrorCode AddGradEqs2D(Mat *A, const PetscInt P_S_inter, const PetscInt n, PetscInt **GI_V_elem,
  PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  MatSetValues(*A, 1 , &GI_V_elem[0][n], P_S_inter, GI_V_elem[2], Gek_x[n], ADD_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[1][n], P_S_inter, GI_V_elem[2], Gek_y[n], ADD_VALUES);  
  return 0;
}

PetscErrorCode AddGradEqs3D(Mat *A, const PetscInt P_S_inter, const PetscInt n, PetscInt **GI_V_elem,
  PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  MatSetValues(*A, 1 , &GI_V_elem[0][n], P_S_inter, GI_V_elem[3], Gek_x[n], ADD_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[1][n], P_S_inter, GI_V_elem[3], Gek_y[n], ADD_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[2][n], P_S_inter, GI_V_elem[3], Gek_z[n], ADD_VALUES);  
  return 0;
}

PetscErrorCode InsertGradEqs2D(Mat *A, const PetscInt P_S_inter, const PetscInt n, PetscInt **GI_V_elem,
  PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  MatSetValues(*A, 1 , &GI_V_elem[0][n], P_S_inter, GI_V_elem[2], Gek_x[n], INSERT_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[1][n], P_S_inter, GI_V_elem[2], Gek_y[n], INSERT_VALUES);  
  return 0;
}

PetscErrorCode InsertGradEqs3D(Mat *A, const PetscInt P_S_inter, const PetscInt n, PetscInt **GI_V_elem,
  PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{
  MatSetValues(*A, 1 , &GI_V_elem[0][n], P_S_inter, GI_V_elem[3], Gek_x[n], INSERT_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[1][n], P_S_inter, GI_V_elem[3], Gek_y[n], INSERT_VALUES);  
  MatSetValues(*A, 1 , &GI_V_elem[2][n], P_S_inter, GI_V_elem[3], Gek_z[n], INSERT_VALUES);  
  return 0;
}


PetscErrorCode ConstViscosity(PetscScalar *Mu, const PetscScalar *mu, StokesVariables_Struct *stokesVariables, 
  const PetscInt *P_S_inter, const PetscInt *inter, const PetscScalar *x, const PetscScalar *y, const PetscScalar *z)
{  *Mu = *mu;   return 0;}

PetscErrorCode CalcViscosityEllipse(PetscScalar *Mu, const PetscScalar *mu, StokesVariables_Struct *stokesVariables, 
  const PetscInt *P_S_inter, const PetscInt *inter, const PetscScalar *x, const PetscScalar *y, const PetscScalar *z)
{
  const PetscScalar mu2 = stokesVariables->mu2;
  const PetscScalar A = stokesVariables->A;
  const PetscScalar B = stokesVariables->B;
  const PetscScalar C = stokesVariables->C;

  *Mu = 0;

  for (PetscInt i=0; i<*P_S_inter;++i) 
  { 
    const PetscInt n = inter[i];
    if (CalcEllipseEquation(x[n], y[n], z[n], A, B, C) >= 0.99)   
      *Mu += mu2;
    else
      *Mu += *mu;
  }
  *Mu = *Mu/(*P_S_inter);

  return 0;
}

PetscErrorCode CalcViscosity2Ellipses(PetscScalar *Mu, const PetscScalar *mu, StokesVariables_Struct *stokesVariables, 
  const PetscInt *P_S_inter, const PetscInt *inter, const PetscScalar *x, const PetscScalar *y, const PetscScalar *z)
{
  const PetscScalar mu2 = stokesVariables->mu2;
  const PetscScalar A  = stokesVariables->A;
  const PetscScalar Bd = stokesVariables->Bd;
  const PetscScalar Bv = stokesVariables->Bv;
  const PetscScalar C  = stokesVariables->C;

  *Mu = 0;

  for (PetscInt i=0; i<*P_S_inter; ++i) 
  { 
    const PetscInt n = inter[i];
    if (y[n] > 0)
    { 
      if (CalcEllipseEquation(x[n], y[n], z[n], A, Bd, C) >= 0.99) 
        *Mu += mu2;
      else
        *Mu += *mu;
    }
    else
    { 
      if (CalcEllipseEquation(x[n], y[n], z[n], A, Bv, C) >= 0.99) 
        *Mu += mu2;
      else
        *Mu += *mu;
    }
  }
  *Mu = *Mu/(*P_S_inter);

  return 0;
}

typedef PetscErrorCode (*AddEqsMom)(Mat *A, Vec *b, const PetscInt U_S_inter, const PetscInt P_S_inter, const PetscInt n, 
  PetscInt **GI_V_elem, PetscScalar **Dek, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z, PetscScalar *Fek_x, 
  PetscScalar *Fek_y, PetscScalar *Fek_z);

typedef PetscErrorCode (*AddEqsCont)(Mat *A, const PetscInt U_S_inter, const PetscInt P_S_inter, const PetscInt n, 
  PetscInt **GI_V_elem, PetscScalar **Mpek, PetscScalar **Gek_xT, PetscScalar **Gek_yT, PetscScalar **Gek_zT);

typedef PetscErrorCode (*AddEqsGrad)(Mat *A, const PetscInt P_S_inter, const PetscInt n, PetscInt **GI_V_elem,
  PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z);

typedef PetscErrorCode (*InsertEqsGrad)(Mat *A, const PetscInt P_S_inter, const PetscInt n, PetscInt **GI_V_elem,
  PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z);

typedef PetscErrorCode (*CalcVisc)(PetscScalar *Mu, const PetscScalar *mu, StokesVariables_Struct *stokesVariables, 
  const PetscInt *P_S_inter, const PetscInt *inter, const PetscScalar *x, const PetscScalar *y, const PetscScalar *z);

AddEqsMom DefineAddMomentumFunction(const PetscInt Dim)
{
	if(Dim == 2)
		return AddMomEqs2D;
	else
		return AddMomEqs3D;
}

AddEqsCont DefineAddContinuityFunction(const PetscInt Dim, const PetscInt Type) 
{
  if (Type == 0)        
    if(Dim == 2)
      return AddContEqs2DNoPress;
    else
      return AddContEqs3DNoPress;
  else if (Type == 1)   
    if(Dim == 2)
      return AddContEqs2DZeroDiagPress;
    else
      return AddContEqs3DZeroDiagPress;
  else                  
    if(Dim == 2)
      return AddContEqs2DMassPress;
    else
      return AddContEqs3DMassPress;
}

AddEqsGrad DefineAddGradientFunction(const PetscInt Dim)
{
  if(Dim == 2)
    return AddGradEqs2D;
  else
    return AddGradEqs3D;
}

InsertEqsGrad DefineInsertGradientFunction(const PetscInt Dim)
{
  if(Dim == 2)
    return InsertGradEqs2D;
  else
    return InsertGradEqs3D;
}

CalcVisc DefineCalcViscosityFunction(const FLUID_VISCOSITY_METHOD ViscMethod)
{
  if(ViscMethod == ONE_CONSTANT)
    return ConstViscosity;
  else if(ViscMethod == TWO_CONSTANTS_ONE_ELLIPSOID)
    return CalcViscosityEllipse;
  else 
    return CalcViscosity2Ellipses;
}


#endif