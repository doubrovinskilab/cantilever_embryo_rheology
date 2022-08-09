
typedef struct {
  FLUID_ELEMENT_TYPE P_elem_type;  
  FLUID_ELEMENT_TYPE U_elem_type;  
  PetscInt     N_inter_glo;  
  PetscInt     N_bound_glo;  
  PetscInt     N_inter_loc;  
  PetscInt     N_bound_loc;  
  PetscInt     N_ghost_loc;    
  PetscInt     P_S_inter;    
  PetscInt     P_S_bound;    
  PetscInt     U_S_inter;    
  PetscInt     U_S_bound;    
  PetscInt     N_Gr;         
  PetscInt     **inter;   
  PetscInt     **bound;   
  PetscInt      *Gr_bound;
  } Elements_Struct;

PetscErrorCode ElementsUpdateInfo(const PetscInt Dim, const FLUID_ELEMENT_TYPE ET, PetscInt *S_inter, PetscInt *S_bound)
{
  if (Dim == 2)
  {
    if ( ET == P1)       *S_inter=3, *S_bound=2;
    else if ( ET == P1b) *S_inter=4, *S_bound=2;
    else if ( ET == P2 ) *S_inter=6, *S_bound=3;
    else if ( ET == Q1 ) *S_inter=4, *S_bound=2;
    else if ( ET == Q2 ) *S_inter=8, *S_bound=3;
    else{
      PetscErrorPrintf("Undefined element type %s in Class_Elements.\n", FLUID_ELEMENT_TYPE_CHAR[ET]);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type.");}
  }
  else if (Dim==3) 
  {
    if ( ET == P1)       *S_inter=4,  *S_bound=3;  
    else if ( ET == P2)  *S_inter=10, *S_bound=6;
    else if ( ET == Q1)  *S_inter=6,  *S_bound=4;
    else if ( ET == P1b) *S_inter=5,  *S_bound=3;  
    else{
      PetscErrorPrintf("Undefined element type %s in Class_Elements.\n", FLUID_ELEMENT_TYPE_CHAR[ET]);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Undefined element type.");}
  }
  else
  {
    PetscErrorPrintf("Incorrect dimension %d in Class_Elements.\n", Dim);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Incorrect dimension.");
  }
  
  return 0; 
}


void ElementsInitialize(Elements_Struct *elements)
{
  elements->P_elem_type = NO_ET;
  elements->U_elem_type = NO_ET;
  elements->N_inter_glo = -1;
  elements->N_bound_glo = -1;
  elements->N_inter_loc = -1;
  elements->N_bound_loc = -1;
  elements->P_S_inter = -1;
  elements->P_S_bound = -1;
  elements->U_S_inter = -1;
  elements->U_S_bound = -1;
   elements->N_Gr = -1;  
}


PetscErrorCode ElementsDefineType(System_Struct *system, Elements_Struct *elements)
{
  if (system->Dim != 2 && system->Dim != 3) {
    PetscErrorPrintf("Error elements dim is not 2 nor 3 (undefined).\n");
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error undefined element dimension.");  }
  if (elements->P_elem_type == NO_ET) {
    PetscErrorPrintf("Error P element type is undefined.\n");
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error undefined element type.");  }
  if (elements->U_elem_type == NO_ET) {
    PetscErrorPrintf("Error U element type is undefined.\n");
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error undefined element type.");  }

  ElementsUpdateInfo(system->Dim, elements->P_elem_type, &elements->P_S_inter, &elements->P_S_bound);
  ElementsUpdateInfo(system->Dim, elements->U_elem_type, &elements->U_S_inter, &elements->U_S_bound);

  return 0;
}     

PetscErrorCode ElementsMallocArrays(System_Struct *system, Elements_Struct *elements)
{
  if(elements->N_inter_loc <= 0)  {
    PetscErrorPrintf("Error N_inter_loc = %d in Elements_Class object rank %d can not be <= 0.\n", elements->N_inter_loc, system->rank);
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error N_inter_loc.");   }
  if(elements->N_bound_loc < 0)   {
    PetscErrorPrintf("Error N_bound_loc = %d in Elements_Class object rank %d can not be <= 0.\n", elements->N_bound_loc, system->rank);
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error N_bound_loc.");   }
  if(elements->U_S_bound <= 0)    {
    PetscErrorPrintf("Error U_S_bound = %d in Elements_Class object rank %d can not be <= 0.\n", elements->U_S_bound, system->rank);
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error U_S_bound.");     }
  if(elements->U_S_inter <= 0)    {
    PetscErrorPrintf("Error U_S_inter = %d in Elements_Class object rank %d can not be <= 0.\n", elements->U_S_inter, system->rank);
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error U_S_inter.");     }

  PetscMalloc1(elements->N_inter_loc, &elements->inter);
  for(PetscInt i = 0; i < elements->N_inter_loc; ++i)   PetscMalloc1(elements->U_S_inter, &elements->inter[i]);
  PetscMalloc1(elements->N_bound_loc, &elements->bound);
  for(PetscInt i = 0; i < elements->N_bound_loc; ++i)   PetscMalloc1(elements->U_S_bound, &elements->bound[i]);
  PetscMalloc1(elements->N_bound_loc, &elements->Gr_bound);
  return 0; 
}

PetscErrorCode ElementsFree(Elements_Struct *elements)
{
  PetscPrintf(PETSC_COMM_WORLD,"  Elements memory is being freed\n");
  for(PetscInt i = 0; i < elements->N_inter_loc; ++i) PetscFree(elements->inter[i]);
  PetscFree(elements->inter);
  for(PetscInt i = 0; i < elements->N_bound_loc; ++i) PetscFree(elements->bound[i]);
  PetscFree(elements->bound);
  PetscFree(elements->Gr_bound);
  return 0;
}

PetscErrorCode ElementsPrintInfo(System_Struct *system, Elements_Struct *elements, PetscBool boolElements) 
{  
  if (boolElements)
  {
    PetscPrintf(PETSC_COMM_WORLD,"  Elements Info:\n");
    PetscPrintf(PETSC_COMM_WORLD,"  Total global internal elements: %d\n", elements->N_inter_glo);
    PetscPrintf(PETSC_COMM_WORLD,"  Total global boundary elements: %d\n", elements->N_bound_glo);

    PetscPrintf(PETSC_COMM_WORLD,"  Internal U elements size: %d\n", elements->U_S_inter);
    PetscPrintf(PETSC_COMM_WORLD,"  Boundary U elements size: %d\n", elements->U_S_bound);

    PetscPrintf(PETSC_COMM_WORLD,"  Internal P elements size: %d\n", elements->P_S_inter);
    PetscPrintf(PETSC_COMM_WORLD,"  Boundary P elements size: %d\n", elements->P_S_bound);

    PetscPrintf(PETSC_COMM_WORLD,"  Total local internal elements:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: %d\n", system->rank, elements->N_inter_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"  Total local boundary elements:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: %d\n", system->rank, elements->N_bound_loc);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
  }
  return 0; 
}