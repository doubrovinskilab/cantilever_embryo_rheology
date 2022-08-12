#ifndef ElementalFunctions_C
#define ElementalFunctions_C

#include "ElementalDelta.c"
#include "ElementalBases.c"
#include "ElementalBasesGradient.c"
#include "ElementalMatrixDiffusion.c"
#include "ElementalMatrixGradient.c"
#include "ElementalMatrixMass.c"
#include "ElementalVectorSource.c"
#include "ElementalVectorBoundaryNeumann.c"
#include "ElementalMatrixBoundaryRobin.c"
#include "ElementalMatrixBoundaryPressure.c"
#include "ElementalMatrixBoundaryPressureNeumann.c"
#include "ElementalMatrixScale.c"
#include "ElementalMatrixScale.c"
#include "ElementalVectorBoundaryNormalTangent.c"


// Defining a function pointers:
typedef PetscErrorCode (*DeltaFunction)( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, PetscScalar *Delta); 
typedef PetscErrorCode (*BaseFunction)( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *phi); 
typedef PetscErrorCode (*GradientBaseFunction)( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, 
  PetscInt **_inter, const PetscScalar xp, const PetscScalar yp, const PetscScalar zp,
  PetscScalar *grad_phi_x, PetscScalar *grad_phi_y, PetscScalar *grad_phi_z); 
typedef PetscErrorCode (*DiffusionFunction)(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
  PetscInt **e_inter, PetscScalar **Dek); 
typedef PetscErrorCode (*GradientFunction)(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, 
  PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z); 
typedef PetscErrorCode (*SourceFunction)(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
 PetscInt **e_inter, PetscScalar *S, PetscScalar *Fek, const PetscBool TypeS);
typedef PetscErrorCode (*MassFunction)(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
  PetscInt **e_inter, PetscScalar **Mek); 
typedef PetscErrorCode (*ScaleFunction)(PetscScalar **Sek, const PetscScalar S); 
typedef PetscErrorCode (*BoundaryNeumannFunction)(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
  PetscInt **e_bound, PetscScalar *Unlp); 
typedef PetscErrorCode (*BoundaryRobinFunction)(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  
  PetscInt **e_bound, PetscScalar **Urlp); 
typedef PetscErrorCode (*BoundaryPressureFunction)(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  
  PetscInt **e_bound, const PetscInt v, PetscScalar **Plp, const PetscScalar *DC); 
typedef PetscErrorCode (*BoundaryNeumannPressureFunction)(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, 
  PetscInt **e_inter, const PetscInt N_inter_loc, PetscScalar **Pnlp, const PetscScalar *DC, const PetscMPIInt rank);
typedef PetscErrorCode (*BoundaryNormalTangentFunction)(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_bound, PetscScalar *N, 
  PetscScalar *T1, PetscScalar *T2, const PetscScalar *DC);

// Defining a error functions:
PetscErrorCode GetDimError_Delta( const PetscInt en, 
  PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter, PetscScalar *Delta)
{ PetscErrorPrintf("Incorrect element dimension for Delta Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Delta Function."); 
  return 0; }

PetscErrorCode GetTypeError_Dek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, 
  PetscScalar **Dek)
{ PetscErrorPrintf("Undefined element type for Diffusion Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Diffusion Function."); 
  return 0; }
PetscErrorCode GetDimError_Dek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, 
  PetscScalar **Dek)
{ PetscErrorPrintf("Incorrect element dimension for Diffusion Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Diffusion Function."); 
  return 0; }

PetscErrorCode GetTypeError_Phi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter,
  const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *phi)
{ PetscErrorPrintf("Undefined element type for Base Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Base Function."); 
  return 0; }
PetscErrorCode GetDimError_Phi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter,
  const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *phi)
{ PetscErrorPrintf("Incorrect element dimension for Base Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Base Function."); 
  return 0; }

PetscErrorCode GetTypeError_GradPhi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter,
  const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *grad_phi_x, PetscScalar *grad_phi_y, 
  PetscScalar *grad_phi_z)
{ PetscErrorPrintf("Undefined element type for Gradient Base Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Gradient Base Function."); 
  return 0; }
PetscErrorCode GetDimError_GradPhi( const PetscInt en, PetscScalar *_x, PetscScalar *_y, PetscScalar *_z, PetscInt **_inter,
  const PetscScalar xp, const PetscScalar yp, const PetscScalar zp, PetscScalar *grad_phi_x, PetscScalar *grad_phi_y, 
  PetscScalar *grad_phi_z)
{ PetscErrorPrintf("Incorrect element dimension for Gradient Base Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Gradient Base Function."); 
  return 0; }

PetscErrorCode GetTypeError_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
  PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{ PetscErrorPrintf("Undefined element type for Gradient Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Gradient Function."); 
  return 0; }
PetscErrorCode GetDimError_Gek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
  PetscInt **e_inter, PetscScalar **Gek_x, PetscScalar **Gek_y, PetscScalar **Gek_z)
{ PetscErrorPrintf("Incorrect element dimension for Gradient Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Gradient Function."); 
  return 0; }

PetscErrorCode GetTypeError_Fek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, 
  PetscInt **e_inter, PetscScalar *S, PetscScalar *Fek, const PetscBool TypeS)
{ PetscErrorPrintf("Undefined element type for Source Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Source Function."); 
  return 0; }
PetscErrorCode GetDimError_Fek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, 
  PetscInt **e_inter, PetscScalar *S, PetscScalar *Fek, const PetscBool TypeS)
{ PetscErrorPrintf("Incorrect element dimension for Source Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Source Function."); 
  return 0; }

PetscErrorCode GetTypeError_Mek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, 
  PetscScalar **Mek)
{ PetscErrorPrintf("Undefined element type for Mass Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Mass Function."); 
  return 0; }
PetscErrorCode GetDimError_Mek(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_inter, 
  PetscScalar **Mek)
{ PetscErrorPrintf("Incorrect element dimension for Mass Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Mass Function."); 
  return 0; }

PetscErrorCode GetTypeError_Scale(PetscScalar **Sek, const PetscScalar S)
{ PetscErrorPrintf("Undefined element type for Scale Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Scale Function."); 
  return 0; }
PetscErrorCode GetDimError_Scale(PetscScalar **Sek, const PetscScalar S)
{ PetscErrorPrintf("Incorrect element dimension for Scale Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Scale Function."); 
  return 0; }

PetscErrorCode GetTypeError_Unlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
  PetscInt **e_bound, PetscScalar *Unlp)
{ PetscErrorPrintf("Undefined element type for Neumann/Robin Boundary Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Neumann/Robin boundary."); 
  return 0; }
PetscErrorCode GetDimError_Unlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,
  PetscInt **e_bound, PetscScalar *Unlp)
{ PetscErrorPrintf("Incorrect element dimension for Neumann/Robin Boundary Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Neumann/Robin boundary."); 
  return 0; }

PetscErrorCode GetTypeError_Urlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  
  PetscInt **e_bound, PetscScalar **Urlp)
{ PetscErrorPrintf("Undefined element type for Robin Boundary Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Robin boundary."); 
  return 0; }
PetscErrorCode GetDimError_Urlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  
  PetscInt **e_bound, PetscScalar **Urlp)
{ PetscErrorPrintf("Incorrect element dimension for Robin Boundary Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Robin boundary."); 
  return 0; }
  
PetscErrorCode GetTypeError_Plp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  
  PetscInt **e_bound, const PetscInt v, PetscScalar **Plp, const PetscScalar *DC)
{ PetscErrorPrintf("Undefined element type for Pressure Boundary Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Pressure Boundary."); 
  return 0; }
PetscErrorCode GetDimError_Plp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  
  PetscInt **e_bound, const PetscInt v, PetscScalar **Plp, const PetscScalar *DC)
{ PetscErrorPrintf("Incorrect element dimension for Pressure Boundary Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Pressure Boundary."); 
  return 0; }

PetscErrorCode GetTypeError_Pnlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, 
  PetscInt **e_inter, const PetscInt N_inter_loc, PetscScalar **Pnlp, const PetscScalar *DC, const PetscMPIInt rank)
{ PetscErrorPrintf("Undefined element type for Neumann Pressure Boundary Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type for Neumann Pressure Boundary."); 
  return 0; }
PetscErrorCode GetDimError_Pnlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, 
  PetscInt **e_inter, const PetscInt N_inter_loc, PetscScalar **Pnlp, const PetscScalar *DC, const PetscMPIInt rank)
{ PetscErrorPrintf("Incorrect element dimension for Neumann Pressure Boundary Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Neumann Pressure Boundary."); 
  return 0; }
  
PetscErrorCode GetDimError_NormTang(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z, PetscInt **e_bound, PetscScalar *N, 
  PetscScalar *T1, PetscScalar *T2, const PetscScalar *DC)
{ PetscErrorPrintf("Incorrect element dimension for Normal-Tangent Boundary Function.\n");
  SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect element dimension for Normal-Tangent Boundary Function."); 
  return 0; }


// Defining functions:
DeltaFunction DefineDeltaFunction(const PetscInt Dim)
{
  if      (Dim == 2)   return Calc2D_Delta;
  else if (Dim == 3)   return Calc3D_Delta;
  else                 return GetDimError_Delta;
}

BaseFunction DefineBaseFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET)
{
  if (Dim == 2)
  {
    if      (ET == P1)     return Calc2DP1_Phi;
    else if (ET == P2)     return Calc2DP2_Phi;
    else if (ET == P1b)    return Calc2DP1b_Phi;
    else                   return GetTypeError_Phi;
  }
  else if (Dim ==3)
  {
    if      (ET == P1)     return Calc3DP1_Phi;
    else if (ET == P2)     return Calc3DP2_Phi;
    else if (ET == P1b)    return Calc3DP1b_Phi;
    else                   return GetTypeError_Phi;
  }
  else
    return GetDimError_Phi;
}

GradientBaseFunction DefineGradientBaseFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET)
{
  if (Dim == 2)
  {
    if      (ET == P1)     return Calc2DP1_GradPhi;
    else if (ET == P2)     return Calc2DP2_GradPhi;
    else if (ET == P1b)    return Calc2DP1b_GradPhi;
    else                   return GetTypeError_GradPhi;
  }
  else if (Dim ==3)
  {
    if      (ET == P1)     return Calc3DP1_GradPhi;
    else if (ET == P2)     return Calc3DP2_GradPhi;
    else if (ET == P1b)    return Calc3DP1b_GradPhi;
    else                   return GetTypeError_GradPhi;
  }
  else
    return GetDimError_GradPhi;
}

DiffusionFunction DefineDiffusionFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET)
{
  if (Dim == 2)
  {
    if      (ET == P1)     return Calc2DP1_Dek;
    else if (ET == P2)     return Calc2DP2_Dek;
    else if (ET == P1b)    return Calc2DP1b_Dek;
    else                   return GetTypeError_Dek;
  }
  else if (Dim ==3)
  {
    if      (ET == P1)     return Calc3DP1_Dek;
    else if (ET == P2)     return Calc3DP2_Dek;
    else if (ET == P1b)    return Calc3DP1b_Dek;
    else                   return GetTypeError_Dek;
  }
  else
    return GetDimError_Dek;
}

GradientFunction DefineGradientFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET_1, const FLUID_ELEMENT_TYPE ET_2)
{
  if (Dim == 2)
  {
    if      (ET_1 == P1  && ET_2 == P1)    return Calc2DP1P1_Gek;
    else if (ET_1 == P2  && ET_2 == P2)    return Calc2DP2P2_Gek;
    else if (ET_1 == P2  && ET_2 == P1)    return Calc2DP2P1_Gek;
    else if (ET_1 == P1  && ET_2 == P2)    return Calc2DP1P2_Gek;
    else if (ET_1 == P1b && ET_2 == P1)    return Calc2DP1bP1_Gek;
    else if (ET_1 == P1  && ET_2 == P1b)   return Calc2DP1P1b_Gek;
    else if (ET_1 == P1b && ET_2 == P1b)   return Calc2DP1bP1b_Gek;
    else                                   return GetTypeError_Gek;
  }
  else if (Dim == 3)
  {
    if      (ET_1 == P1  && ET_2 == P1)    return Calc3DP1P1_Gek;
    else if (ET_1 == P2  && ET_2 == P2)    return Calc3DP2P2_Gek;
    else if (ET_1 == P2  && ET_2 == P1)    return Calc3DP2P1_Gek;
    else if (ET_1 == P1  && ET_2 == P2)    return Calc3DP1P2_Gek;
    else if (ET_1 == P1b && ET_2 == P1)    return Calc3DP1bP1_Gek;
    else if (ET_1 == P1  && ET_2 == P1b)   return Calc3DP1P1b_Gek;
    else if (ET_1 == P1b && ET_2 == P1b)   return Calc3DP1bP1b_Gek;
    else                                   return GetTypeError_Gek;
  }
  else
    return GetDimError_Gek;
}

SourceFunction DefineSourceFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET)
{
  if (Dim == 2)
  {
    if      (ET == P1)  return Calc2DP1_Fek;
    else if (ET == P2)  return Calc2DP2_Fek;
    else if (ET == P1b) return Calc2DP1b_Fek;
    else                return GetTypeError_Fek;
  }
  else if (Dim ==3)
  {
    if      (ET == P1)  return Calc3DP1_Fek;
    else if (ET == P2)  return Calc3DP2_Fek;
    else if (ET == P1b) return Calc3DP1b_Fek;
    else                return GetTypeError_Fek;
  }
  else
    return GetDimError_Fek;
}

MassFunction DefineMassFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET)
{
  if (Dim == 2)
  {
    if      (ET == P1)     return Calc2DP1_Mek;
    else if (ET == P2)     return Calc2DP2_Mek;
    else if (ET == P1b)    return Calc2DP1b_Mek;
    else                   return GetTypeError_Mek;
  }
  else if (Dim ==3)
  {
    if      (ET == P1)     return Calc3DP1_Mek;
    else if (ET == P2)     return Calc3DP2_Mek;
    else if (ET == P1b)    return Calc3DP1b_Mek;
    else                   return GetTypeError_Mek;
  }
  else
    return GetDimError_Mek;
}

ScaleFunction DefineScaleFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET_1, const FLUID_ELEMENT_TYPE ET_2)
{
  if (Dim == 2)
  {
    if      (ET_1 == P1  && ET_2 == P1)    return Calc2DP1P1_Scale;
    else if (ET_1 == P2  && ET_2 == P2)    return Calc2DP2P2_Scale;
    else if (ET_1 == P2  && ET_2 == P1)    return Calc2DP2P1_Scale;
    else if (ET_1 == P1  && ET_2 == P2)    return Calc2DP1P2_Scale;
    else if (ET_1 == P1b && ET_2 == P1)    return Calc2DP1bP1_Scale;
    else if (ET_1 == P1  && ET_2 == P1b)   return Calc2DP1P1b_Scale;
    else if (ET_1 == P1b && ET_2 == P1b)   return Calc2DP1bP1b_Scale;
    else                                   return GetTypeError_Scale;
  }
  else if (Dim == 3)
  {
    if      (ET_1 == P1  && ET_2 == P1)    return Calc3DP1P1_Scale;
    else if (ET_1 == P2  && ET_2 == P2)    return Calc3DP2P2_Scale;
    else if (ET_1 == P2  && ET_2 == P1)    return Calc3DP2P1_Scale;
    else if (ET_1 == P1  && ET_2 == P2)    return Calc3DP1P2_Scale;
    else if (ET_1 == P1b && ET_2 == P1)    return Calc3DP1bP1_Scale;
    else if (ET_1 == P1  && ET_2 == P1b)   return Calc3DP1P1b_Scale;
    else if (ET_1 == P1b && ET_2 == P1b)   return Calc3DP1bP1b_Scale;
    else                                   return GetTypeError_Scale;
  }
  else
    return GetDimError_Scale;
}

BoundaryNeumannFunction DefineBoundaryNeumannFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET)
{
  if (Dim == 2)
  {
    if      (ET == P1)  return Calc2DP1_Unlp;
    else if (ET == P2)  return Calc2DP2_Unlp;
    else if (ET == P1b) return Calc2DP1b_Unlp;
    else                return GetTypeError_Unlp;
  }
  else if (Dim ==3)
  {
    if      (ET == P1)  return Calc3DP1_Unlp;
    else if (ET == P2)  return Calc3DP2_Unlp;
    else if (ET == P1b) return Calc3DP1b_Unlp;
    else                return GetTypeError_Unlp;
  }
  else
    return GetDimError_Unlp;
}

BoundaryRobinFunction DefineBoundaryRobinFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET)
{
  if (Dim == 2)
  {
    if      (ET == P2)  return Calc2DP2_Urlp;
    else if (ET == P1b) return Calc2DP1b_Urlp;
    else                return GetTypeError_Urlp;
  }
  else if (Dim ==3)
  {
    if      (ET == P2)  return Calc3DP2_Urlp;
    else if (ET == P1b) return Calc3DP1b_Urlp;
    else                return GetTypeError_Urlp;
  }
  else
    return GetDimError_Urlp;
}

BoundaryPressureFunction DefineBoundaryPressureFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET_1, const FLUID_ELEMENT_TYPE ET_2)
{
  if (Dim == 2)
  {
    if      (ET_1 == P2 && ET_2 == P1)    return Calc2DP2P1_Plp;
    else if (ET_1 == P1b && ET_2 == P1)   return Calc2DP1bP1_Plp;
    else                                  return GetTypeError_Plp;
  }
  else if (Dim ==3)
  {
    if      (ET_1 == P2 && ET_2 == P1)    return Calc3DP2P1_Plp;
    else if (ET_1 == P1b && ET_2 == P1)   return Calc3DP1bP1_Plp;
    else                                  return GetTypeError_Plp;
  }
  else
    return GetDimError_Plp;
}


BoundaryNeumannPressureFunction DefineBoundaryNeumannPressureFunction(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET)
{
  if (Dim == 2)
  {
    if      (ET == P1)  return Calc2DP1_Pnlp;
    else                return GetTypeError_Pnlp;
  }
  else if (Dim ==3)
  {
    if      (ET == P1)  return Calc3DP1_Pnlp;
    else                return GetTypeError_Pnlp;
  }
  else
    return GetDimError_Pnlp;
}

BoundaryNormalTangentFunction DefineBoundaryNormalTangentFunction(const PetscInt Dim)
{
  if (Dim == 2)
    return Calc2D_NormTang;
  else if (Dim ==3)
    return Calc3D_NormTang;
  else
    return GetDimError_NormTang;
}
#endif