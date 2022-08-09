
typedef PetscErrorCode (*SearchForForceFunction)(System_Struct *system, Solid_Struct *solid); 
typedef PetscErrorCode (*SearchForForceFunction2)(Solid_Struct *solid); 

typedef struct 
{
  SearchForForceFunction LE_func; 
  SearchForForceFunction SF_func; 
  SearchForForceFunction2 BC_InsertForce_func;          
  SearchForForceFunction2 BC_AddForce_func;             
  SearchForForceFunction2 BC_AddForceEllipse_func;      
  SearchForForceFunction  BC_SpeedForce_func;           
  SearchForForceFunction  BC_SpeedForceEllipse_func;    
  SearchForForceFunction2 BC_DistFixed_func;              
} ForcesActive_Struct; 


SearchForForceFunction DefineSearchForLinearElasticFunction(Solid_Struct *solid)
{
  if(solid->forces.LE.Type == FORCE_LE_FIXED_DISP) {
    return LinearElasticity;}
  else if(solid->forces.LE.Type == FORCE_LE_INITIAL_DISP){
    return LinearElasticityOrigin; }
  else if(solid->forces.LE.Type == FORCE_LE_INITIAL_DISP_MULTI){
    return LinearElasticityOriginMultiK; }
  else if(solid->forces.LE.Type == FORCE_LE_INITIAL_DISP_COND){
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! FORCE_LE_INITIAL_DISP_COND not defined");  
    return LinearElasticityNull; }
  else {
    return LinearElasticityNull;  }
}


SearchForForceFunction DefineSearchForSoftForceFunction(Solid_Struct *solid)
{
  if (solid->forces.SF.Type == FORCE_SF_ELLIPSOID){
    return SoftForceOnOneEllipsoid;}
  else if (solid->forces.SF.Type == FORCE_SF_TWO_ELLIPSOIDS){
    return SoftForceOnTwoEllipsoids;}
  else {
    return SoftForceNull;}

}


SearchForForceFunction2 DefineSearchForBFInsertForceFunction(Solid_Struct *solid)
{
  PetscBool BC_InsertForce_pts = PETSC_FALSE,  BC_InsertForce_lns = PETSC_FALSE,  BC_InsertForce_trs = PETSC_FALSE;
  if (solid->N_bound != -1 )  
  {
    for (PetscInt bp = 0; bp<solid->N_bound_pts; bp++) {
      const PetscInt g = solid->bound_group_pts[bp];  
      if (solid->Gr_type[g] == FORCE_INSERT)   BC_InsertForce_pts = PETSC_TRUE; }

    for (PetscInt bl = 0; bl<solid->N_bound_lns; bl++) {
      const PetscInt g = solid->bound_group_lns[bl];  
      if (solid->Gr_type[g] == FORCE_INSERT) BC_InsertForce_lns = PETSC_TRUE; }

    for (PetscInt bt = 0; bt<solid->N_bound_trs; bt++) {
      const PetscInt g = solid->bound_group_trs[bt];  
      if (solid->Gr_type[g] == FORCE_INSERT) BC_InsertForce_trs = PETSC_TRUE;}

    if ((PetscInt) (BC_InsertForce_pts+BC_InsertForce_lns+BC_InsertForce_trs)> 1){
      return BoundaryForce_InsertForce;}
    else if (BC_InsertForce_pts){
      return BoundaryForce_InsertForcePoints;}
    else if (BC_InsertForce_lns){
      return BoundaryForce_InsertForceLines;}
    else if (BC_InsertForce_trs){
      return BoundaryForce_InsertForceTriangles;}
    else {
      return BoundaryForce_InsertForceNull;}
  }
  else {
    return BoundaryForce_InsertForceNull;}
}


SearchForForceFunction2 DefineSearchForBFAddForceFunction(Solid_Struct *solid)
{
  PetscBool BC_AddForce_pts = PETSC_FALSE,  BC_AddForce_lns = PETSC_FALSE,  BC_AddForce_trs = PETSC_FALSE;
  if (solid->N_bound != -1 )  
  {
    for (PetscInt bp = 0; bp<solid->N_bound_pts; bp++) {
      const PetscInt g = solid->bound_group_pts[bp];  
      if (solid->Gr_type[g] == FORCE_ADD)   BC_AddForce_pts = PETSC_TRUE; }

    for (PetscInt bl = 0; bl<solid->N_bound_lns; bl++) {
      const PetscInt g = solid->bound_group_lns[bl];  
      if (solid->Gr_type[g] == FORCE_ADD) BC_AddForce_lns = PETSC_TRUE; }

    for (PetscInt bt = 0; bt<solid->N_bound_trs; bt++) {
      const PetscInt g = solid->bound_group_trs[bt];  
      if (solid->Gr_type[g] == FORCE_ADD) BC_AddForce_trs = PETSC_TRUE;}

    if ((PetscInt) (BC_AddForce_pts+BC_AddForce_lns+BC_AddForce_trs)> 1){
      return BoundaryForce_AddForce;}
    else if (BC_AddForce_pts){
      return BoundaryForce_AddForcePoints;}
    else if (BC_AddForce_lns){
      return BoundaryForce_AddForceLines;}
    else if (BC_AddForce_trs){
      return BoundaryForce_AddForceTriangles;}
    else {
      return BoundaryForce_AddForceNull;}
  }
  else {
    return BoundaryForce_AddForceNull;}
}

SearchForForceFunction2 DefineSearchForBFAddForceEllipseFunction(Solid_Struct *solid)
{
  PetscBool BC_AddForceEllipse_pts = PETSC_FALSE,  BC_AddForceEllipse_lns = PETSC_FALSE,  BC_AddForceEllipse_trs = PETSC_FALSE;
  if (solid->N_bound != -1 )  
  {
    for (PetscInt bp = 0; bp<solid->N_bound_pts; bp++) {
      const PetscInt g = solid->bound_group_pts[bp];  
      if (solid->Gr_type[g] == FORCE_ADD_ELLIPSE)   BC_AddForceEllipse_pts = PETSC_TRUE; }

    for (PetscInt bl = 0; bl<solid->N_bound_lns; bl++) {
      const PetscInt g = solid->bound_group_lns[bl];  
      if (solid->Gr_type[g] == FORCE_ADD_ELLIPSE) BC_AddForceEllipse_lns = PETSC_TRUE; }

    for (PetscInt bt = 0; bt<solid->N_bound_trs; bt++) {
      const PetscInt g = solid->bound_group_trs[bt];  
      if (solid->Gr_type[g] == FORCE_ADD_ELLIPSE) BC_AddForceEllipse_trs = PETSC_TRUE;}

    if ((PetscInt) (BC_AddForceEllipse_pts+BC_AddForceEllipse_lns+BC_AddForceEllipse_trs)> 1){    
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error in finding the right BoundaryForce_AddForceEllipse function");
      return BoundaryForce_AddForceEllipsePoints;}
    else if (BC_AddForceEllipse_pts){
      return BoundaryForce_AddForceEllipsePoints;}
    else if (BC_AddForceEllipse_lns){
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error in finding the right BoundaryForce_AddForceEllipse function");
      return BoundaryForce_AddForceEllipseNull;}
    else if (BC_AddForceEllipse_trs){
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error in finding the right BoundaryForce_AddForceEllipse function");
      return BoundaryForce_AddForceEllipseNull;}
    else {
      return BoundaryForce_AddForceEllipseNull;}
  }
  else {
    return BoundaryForce_AddForceEllipseNull;}
}


SearchForForceFunction DefineSearchForBFSpeedForceFunction(Solid_Struct *solid)
{
  if (solid->N_bound != -1 )  
  {
    PetscBool Found = PETSC_FALSE;

    for(PetscInt j=0; j< solid->N_bound_pts; j++) {   
      const PetscInt g = solid->bound_group_pts[j];   
      if (solid->Gr_type[g] == FORCE_SPEED) Found = PETSC_TRUE;}

    for(PetscInt j=0; j< solid->N_bound_lns; j++) {   
      const PetscInt g = solid->bound_group_lns[j];   
      if (solid->Gr_type[g] == FORCE_SPEED) Found = PETSC_TRUE; }

    for(PetscInt j=0; j< solid->N_bound_trs; j++) {   
      const PetscInt g = solid->bound_group_trs[j];   
      if (solid->Gr_type[g] == FORCE_SPEED) Found = PETSC_TRUE; }

    if (Found){
      return BoundaryForce_SpeedForce;}
    else {
      return BoundaryForce_SpeedForceNull;}
  }
  else {
    return BoundaryForce_SpeedForceNull;}
}

SearchForForceFunction DefineSearchForBFSpeedForceEllipseFunction(Solid_Struct *solid)
{
  if (solid->N_bound != -1 )  
  {    
    PetscBool Found = PETSC_FALSE;

    for(PetscInt j=0; j< solid->N_bound_pts; j++) {   
      const PetscInt g = solid->bound_group_pts[j];   
      if (solid->Gr_type[g] == FORCE_SPEED_ELLIPSE)  Found = PETSC_TRUE;}

    for(PetscInt j=0; j< solid->N_bound_lns; j++) {   
      const PetscInt g = solid->bound_group_lns[j];   
      if (solid->Gr_type[g] == FORCE_SPEED_ELLIPSE)  Found = PETSC_TRUE; }

    for(PetscInt j=0; j< solid->N_bound_trs; j++) {   
      const PetscInt g = solid->bound_group_trs[j];   
      if (solid->Gr_type[g] == FORCE_SPEED_ELLIPSE)  Found = PETSC_TRUE; }

    if (Found){
      return BoundaryForce_SpeedForceEllipse;}
    else {
      return BoundaryForce_SpeedForceEllipseNull;}
  }
  else {
    return BoundaryForce_SpeedForceEllipseNull;}
}


SearchForForceFunction2 DefineSearchForBFDistFixedFunction(Solid_Struct *solid)
{
  if (solid->N_bound != -1 )  
  {    
    PetscBool Found = PETSC_FALSE;

    for (PetscInt bp = 0; bp<solid->N_bound_pts; bp++) {
      const PetscInt g = solid->bound_group_pts[bp];  
      if (solid->Gr_type[g] == DIST_FIXED) Found = PETSC_TRUE; }

    for (PetscInt bl = 0; bl<solid->N_bound_lns; bl++) {
      const PetscInt g = solid->bound_group_lns[bl];  
      if (solid->Gr_type[g] == DIST_FIXED) Found = PETSC_TRUE; }

    for (PetscInt bt = 0; bt<solid->N_bound_trs; bt++) {
      const PetscInt g = solid->bound_group_trs[bt];  
      if (solid->Gr_type[g] == DIST_FIXED) Found = PETSC_TRUE; }

    if (Found){
      return BoundaryForce_DistFixed;}
    else {
      return BoundaryForce_DistFixedNull;}
  } 
  else {
    return BoundaryForce_DistFixedNull;}
}



PetscErrorCode ForcesInitialize(System_Struct *system, Solid_Struct *solid, ForcesActive_Struct *active_forces)
{
  ForcesMallocPointArrays(solid);
  active_forces->LE_func = DefineSearchForLinearElasticFunction(solid);
  LinearElasticityInitialize(solid);
  active_forces->SF_func = DefineSearchForSoftForceFunction(solid);
  SoftForceInitialization(system, solid);
  
  active_forces->BC_InsertForce_func = DefineSearchForBFInsertForceFunction(solid);
  active_forces->BC_AddForce_func = DefineSearchForBFAddForceFunction(solid);
  active_forces->BC_AddForceEllipse_func = DefineSearchForBFAddForceEllipseFunction(solid);
  active_forces->BC_SpeedForce_func = DefineSearchForBFSpeedForceFunction(solid);
  BoundaryForce_SpeedForceInitialize(solid);
  active_forces->BC_SpeedForceEllipse_func = DefineSearchForBFSpeedForceEllipseFunction(solid);
  BoundaryForce_SpeedForceEllipseInitialize(solid);
  active_forces->BC_DistFixed_func = DefineSearchForBFDistFixedFunction(solid);
  return 0;
}


PetscErrorCode ForcesCalculate(System_Struct *system, Solid_Struct *solid, ForcesActive_Struct *active_forces)
{
  ForcesReset(solid);
  active_forces->LE_func(system, solid);
  active_forces->SF_func(system, solid);
   
  active_forces->BC_InsertForce_func(solid);
  active_forces->BC_AddForce_func(solid);
  active_forces->BC_AddForceEllipse_func(solid);
  active_forces->BC_SpeedForce_func(system, solid);
  active_forces->BC_SpeedForceEllipse_func(system, solid);
  active_forces->BC_DistFixed_func(solid);
  return 0;
}
