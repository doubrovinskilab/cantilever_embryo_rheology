
PetscPrintf(PETSC_COMM_WORLD,"Reading parameters from input file:\n"); 
PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------\n"); 

// System Variables
//=================================================================================================
PetscPrintf(PETSC_COMM_WORLD,"  Reading system parameters:\n"); 
SearchInputAsciiFileForEnumSystemType(system.SYSInputFile,     "System_Type",             &system.SYSType);  
SearchInputAsciiFileForStr(system.SYSInputFile,                "System_OutputFolder",      system.SYSOutputFolder, PETSC_FALSE);
SearchInputAsciiFileForInt(system.SYSInputFile,                "System_Dim",              &system.Dim);   
if (system.SYSType == IBM || system.SYSType == LINEAR_ELASTICITY)
{
  SearchInputAsciiFileForScalar(system.SYSInputFile,           "System_end_time",         &system.end_time);               
  SearchInputAsciiFileForScalar(system.SYSInputFile,           "System_dtime_save",       &system.dtime_save);   // Interval to save time
}

if (system.SYSType == LINEAR_ELASTICITY)
  SearchInputAsciiFileForScalar(system.SYSInputFile,  "Solid_dt",  &system.SD_dt);

// Fluid Variables
//=================================================================================================
if (system.SYSType != LINEAR_ELASTICITY)
{
  // Fluid Solver & Mesh Options
  PetscPrintf(PETSC_COMM_WORLD,"  Reading fluid parameters:\n"); 
  SearchInputAsciiFileForEnumFLSolverType(system.SYSInputFile,   "Fluid_Solver_Type",     &system.FLSolver_Type);  
  SearchInputAsciiFileForEnumFLSolverMethod(system.SYSInputFile, "Fluid_Solver_Method",   &system.FLSolver_Method);  
  SearchInputAsciiFileForEnumMeshType(system.SYSInputFile,       "Fluid_Mesh_Type",       &system.FLMesh_Type);
  SearchInputAsciiFileForStr(system.SYSInputFile,                "Fluid_Mesh_File",        system.FLMesh_File, PETSC_FALSE);

  if (system.FLSolver_Type == EXPLICIT_LUMPED || system.FLSolver_Type == IMPLICIT_EXPLICIT_LUMPED)
  {
    SearchInputAsciiFileForInt(system.SYSInputFile,              "Fluid_Exp_N",    &system.FLExp_N);
    SearchInputAsciiFileForScalar(system.SYSInputFile,           "Fluid_Exp_dt",   &system.FLExp_dt);
    SearchInputAsciiFileForScalar(system.SYSInputFile,           "Fluid_Exp_eps",  &system.FLExp_eps);
    SearchInputAsciiFileForScalar(system.SYSInputFile,           "Fluid_Exp_rho",  &system.FLExp_rho);
  }

  if (system.FLSolver_Type == IMPLICIT || system.FLSolver_Type == IMPLICIT_EXPLICIT_LUMPED)
  {
    SearchInputAsciiFileForScalar(system.SYSInputFile,           "Fluid_Imp_dt",   &system.FLImp_dt);
    SearchInputAsciiFileForScalar(system.SYSInputFile,           "Fluid_Imp_eps",  &system.FLImp_eps);
  }

  if (system.FLSolver_Type == IMPLICIT_EXPLICIT_LUMPED)
    SearchInputAsciiFileForScalar(system.SYSInputFile,           "Fluid_ImpExp_rTol",  &system.FLImpExp_rTol);

  // Fluid Variable & Element Option
  SearchInputAsciiFileForEnumFLElemType(system.SYSInputFile,     "Fluid_ElemClass_P",      &elements.P_elem_type);
  SearchInputAsciiFileForEnumFLElemType(system.SYSInputFile,     "Fluid_ElemClass_U",      &elements.U_elem_type);

  SearchInputAsciiFileForScalar(system.SYSInputFile,             "Fluid_Variables_mu",     &stokesVariables.mu);
  SearchInputAsciiFileForEnumFLViscMethod(system.SYSInputFile,   "Fluid_Calc_Viscosity",   &stokesVariables.ViscMethod);
  if (stokesVariables.ViscMethod != ONE_CONSTANT) 
  {
    SearchInputAsciiFileForScalar(system.SYSInputFile,           "Fluid_Variables_mu2",    &stokesVariables.mu2);
    SearchInputAsciiFileForScalar(system.SYSInputFile,           "Fluid_Ellipsoid_A",      &stokesVariables.A);
    if (stokesVariables.ViscMethod == TWO_CONSTANTS_ONE_ELLIPSOID){
      SearchInputAsciiFileForScalar(system.SYSInputFile,         "Fluid_Ellipsoid_B",      &stokesVariables.B); }
    else if (stokesVariables.ViscMethod == TWO_CONSTANTS_TWO_ELLIPSOIDS){
      SearchInputAsciiFileForScalar(system.SYSInputFile,         "Fluid_Ellipsoid_Bd",      &stokesVariables.Bd);
      SearchInputAsciiFileForScalar(system.SYSInputFile,         "Fluid_Ellipsoid_Bv",      &stokesVariables.Bv);
    }
    if (system.Dim == 3)  SearchInputAsciiFileForScalar(system.SYSInputFile,   "Fluid_Ellipsoid_C",   &stokesVariables.C);
    else                  stokesVariables.C = 1.0;
  }
  SearchInputAsciiFileFor3Scalar(system.SYSInputFile,  "Fluid_Domain_Center",  nodes.DC);
}

// Solid Variables
//=================================================================================================
if (system.SYSType == IBM || system.SYSType == LINEAR_ELASTICITY)
{
  PetscPrintf(PETSC_COMM_WORLD,"  Reading solid parameters:\n"); 
  SearchInputAsciiFileForStr(system.SYSInputFile,              "Solid_Mesh_File",  system.SDMesh_File, PETSC_FALSE);
  SearchInputAsciiFileForEnumMeshType(system.SYSInputFile,     "Solid_Mesh_Type",  &system.SDMesh_Type);
  SearchInputAsciiFileForBoolIfFound(system.SYSInputFile,      "Solid_Mesh_ConvTrs2Lns",    &system.SDMesh_ConvTrs2Lns);

  // Linear Elastic Forces
  { 
    SearchInputAsciiFileForEnumSDForceLETypeIfFound(system.SYSInputFile,"Solid_Force_LE_Type",    &solid.forces.LE.Type);
    if (solid.forces.LE.Type != FORCE_LE_NULL)
      SearchInputAsciiFileForScalar(system.SYSInputFile,         "Solid_Force_LE_K",  &solid.forces.LE.K);
    if (solid.forces.LE.Type == FORCE_LE_INITIAL_DISP_MULTI)
      SearchInputAsciiFileForScalarDataFile(system.SYSInputFile, "Solid_Force_LE_Multi_File", &system, &solid.forces.LE.M_size, &solid.forces.LE.M);
    if (solid.forces.LE.Type == FORCE_LE_FIXED_DISP)
      SearchInputAsciiFileForScalar(system.SYSInputFile,         "Solid_Force_LE_D0",   &solid.forces.LE.L0);
  }


  // Soft Forces
  {
    SearchInputAsciiFileForEnumSDForceSFTypeIfFound(system.SYSInputFile, "Solid_Force_SF_Type",    &solid.forces.SF.Type);
    if (solid.forces.SF.Type == FORCE_SF_ELLIPSOID || solid.forces.SF.Type == FORCE_SF_TWO_ELLIPSOIDS || solid.forces.SF.Type == FORCE_SF_ELLIPSOID_REMOVE_NORMAL_FORCE)
    {
      SearchInputAsciiFileForScalar(system.SYSInputFile,   "Solid_Force_SF_Ellipsoid_A",   &solid.forces.SF.A);
      if (solid.forces.SF.Type == FORCE_SF_ELLIPSOID || solid.forces.SF.Type == FORCE_SF_ELLIPSOID_REMOVE_NORMAL_FORCE)
        SearchInputAsciiFileForScalar(system.SYSInputFile, "Solid_Force_SF_Ellipsoid_B",   &solid.forces.SF.B);
      else if (solid.forces.SF.Type == FORCE_SF_TWO_ELLIPSOIDS)  
      {
        SearchInputAsciiFileForScalar(system.SYSInputFile, "Solid_Force_SF_Ellipsoid_Bd",  &solid.forces.SF.Bd);
        SearchInputAsciiFileForScalar(system.SYSInputFile, "Solid_Force_SF_Ellipsoid_Bv",  &solid.forces.SF.Bv);
      }
      if (system.Dim == 3)  SearchInputAsciiFileForScalar(system.SYSInputFile,  "Solid_Force_SF_Ellipsoid_C",    &solid.forces.SF.C);
      else                  solid.forces.SF.C = 1.0;

      SearchInputAsciiFileForScalar(system.SYSInputFile,   "Solid_Force_SF_Ellipsoid_ds",   &solid.forces.SF.ds);
      SearchInputAsciiFileForScalar(system.SYSInputFile,   "Solid_Force_SF_Ellipsoid_sf",   &solid.forces.SF.sf);
    }
  }

  // Constrain
  {
    SearchInputAsciiFileForEnumSDConstrainTypeIfFound(system.SYSInputFile,"Solid_Constrain_Type",  &solid.constrain.Type);
    if (  solid.constrain.Type == CONSTRAIN_ELLIPSOID ||  solid.constrain.Type == CONSTRAIN_TWO_ELLIPSOIDS)
    {
      SearchInputAsciiFileForScalar(system.SYSInputFile,    "Solid_Constrain_Ellipsoid_A",    &solid.constrain.A);
      if ( solid.constrain.Type == CONSTRAIN_ELLIPSOID)
        SearchInputAsciiFileForScalar(system.SYSInputFile,  "Solid_Constrain_Ellipsoid_B",    &solid.constrain.B);
      else if (solid.constrain.Type == CONSTRAIN_TWO_ELLIPSOIDS)  {
        SearchInputAsciiFileForScalar(system.SYSInputFile,  "Solid_Constrain_Ellipsoid_Bd",   &solid.constrain.Bd);
        SearchInputAsciiFileForScalar(system.SYSInputFile,  "Solid_Constrain_Ellipsoid_Bv",   &solid.constrain.Bv); }

      if (system.Dim == 3)  SearchInputAsciiFileForScalar(system.SYSInputFile, "Solid_Constrain_Ellipsoid_C",    &solid.constrain.C);
      else                  solid.constrain.C = 1.0;
    }
    else if (solid.constrain.Type == CONSTRAIN_AXIS)
    { 
      SearchInputAsciiFileForInt(system.SYSInputFile,     "Solid_Constrain_Axis_Index",   &solid.constrain.axis_i);
      SearchInputAsciiFileForScalar(system.SYSInputFile,  "Solid_Constrain_Axis_Value",   &solid.constrain.axis_v);
    }
    if (solid.constrain.Type != CONSTRAIN_NULL)
      SearchInputAsciiFileForScalar(system.SYSInputFile,  "Solid_Constrain_ds",      &solid.constrain.ds);
  }

  // Make Solid Velocity Divergent Free
  SearchInputAsciiFileForBoolIfFound(system.SYSInputFile,   "Solid_Div_Free",    &system.SDDivFree_Bool);
  SearchInputAsciiFileForIntIfFound(system.SYSInputFile,    "Solid_Dim",         &solid.Dim);   
  if (solid.Dim == -1) solid.Dim = system.Dim;

  if (system.SYSType == LINEAR_ELASTICITY)
    SearchInputAsciiFileForScalar(system.SYSInputFile,     "Solid_Damping_C",   &solid.forces.c);

  // Reload Data
  {
    SearchInputAsciiFileForBool(system.SYSInputFile,        "Solid_Reload",         &system.SDReload);
    if(system.SDReload){
      SearchInputAsciiFileForStr(system.SYSInputFile,       "Solid_Reload_File",          system.SDReload_File, PETSC_FALSE);
      SearchInputAsciiFileForScalar(system.SYSInputFile,    "Solid_Reload_Time",         &system.SDReload_time);
      SearchInputAsciiFileForInt(system.SYSInputFile,       "Solid_Reload_Saved_Step",   &system.SDReload_saved_step);
    }
  }
}

// Check for Incompatibility of the options:
//==============================================
if (system.Dim != 2 && system.Dim != 3)
{ PetscErrorPrintf("Incorrect dimenstion '%d'.\n", system.Dim);
  SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Incorrect dimension.");  }

if (system.SYSType != LINEAR_ELASTICITY)
  if (system.FLSolver_Type == EXPLICIT_LUMPED && (system.FLSolver_Method == DIRECT || system.FLSolver_Method == NEST || system.FLSolver_Method == MONO))
  { PetscErrorPrintf("%s can not work with %s method.\n", FLUID_SOLVER_TYPE_CHAR[system.FLSolver_Type], FLUID_SOLVER_METHOD_CHAR[system.FLSolver_Method]);
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Incorrect method.");  }

if (system.rank == 0) PetscMkdir(system.SYSOutputFolder);

PetscPrintf(PETSC_COMM_WORLD,"\n"); 