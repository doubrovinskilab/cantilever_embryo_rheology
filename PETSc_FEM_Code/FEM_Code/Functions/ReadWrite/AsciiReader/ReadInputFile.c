#ifndef ReadInputFile
#define ReadInputFile


PetscErrorCode FindSizeOfDataFile(const char* file_name, System_Struct *system)
{
  PetscErrorCode ierr;
  FILE *fd;
  char line[256]="";
  PetscInt array_size;
  size_t len;

  if (system->rank == 0)
  {
    fd   = fopen(file_name,"r");

    if(fd == NULL) 
    { PetscErrorPrintf("File %s is not found.\n", file_name);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

    array_size = 0;
    
    while (PETSC_TRUE)
    {
      ReadFileLine(fd, 256,line); 
      PetscStrlen(line, &len);
      if ( len==0)        break; 
      if ( line[0] == 0)  break; 
      array_size++;
    }
    
    ierr = fclose(fd);
    if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");
  }

  if (system->size>1) {
    ierr = MPI_Bcast(&array_size, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); PetscBarrier(NULL); }

  return array_size;
}

PetscErrorCode ReadIntDataFromFile(const char* file_name, const PetscInt array_size, PetscInt* array_int, System_Struct *system)
{
  PetscErrorCode ierr;
  FILE *fd;
  char line[256]="";
  char data2[256];

  if (system->rank == 0)
  {
    fd   = fopen(file_name,"r");

    if(fd == NULL) 
    { PetscErrorPrintf("File %s is not found.\n", file_name);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

    
    for(PetscInt i=0; i < array_size; i++) 
    {
      
      ReadFileLine(fd, 256,line); 
      RemoveNewLine(line, data2); 
      
      ierr = PetscOptionsStringToInt(data2, &array_int[i]);      CHKERRQ(ierr);
    }
    
    ierr = fclose(fd);
    if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");
  }

  if (system->size>1) {
    ierr = MPI_Bcast(array_int, array_size, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); PetscBarrier(NULL); }

  return 0;
}
PetscErrorCode ReadScalarDataFromFile(const char* file_name, const PetscInt array_size, PetscScalar* array_scalar, System_Struct *system)
{
  PetscErrorCode ierr;
  FILE *fd;
  char line[256]="";
  char data2[256];

  if (system->rank == 0)
  {
    fd   = fopen(file_name,"r");

    if(fd == NULL) 
    { PetscErrorPrintf("File %s is not found.\n", file_name);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

    
    for(PetscInt i=0; i < array_size; i++) 
    {
      ReadFileLine(fd, 256,line); 
      RemoveNewLine(line, data2); 
      
      ierr = PetscOptionsStringToScalar(data2, &array_scalar[i]);      CHKERRQ(ierr); 
    }
    
    ierr = fclose(fd);
    if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");
  }

  if (system->size>1) {
    ierr = MPI_Bcast(array_scalar, array_size, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); PetscBarrier(NULL); }

  return 0;
}

PetscErrorCode SearchInputAsciiFile(const char* file_name, const char str_search[], char str_value[],const PetscBool BoolLower)
{
  /* Searches for string in a file. If BoolLower is true ==> change to lower case. */
  PetscErrorCode ierr = 0;
  PetscInt        i;
  PetscBool   found = PETSC_FALSE;
  char data[256]="";
  char **data2;
  char **data3;
  PetscInt j, k;
  char str_value_space[256]="";
  char *str_value_temp;

  FILE *file;
  ierr = PetscFOpen(PETSC_COMM_WORLD, file_name,"r", &file); CHKERRQ(ierr);
  for (i=0; i < 1000; i++)
  {
    PetscSynchronizedFGets(PETSC_COMM_WORLD, file, 256, data);
    PetscStrInList(str_search, data,' ', &found);
    if (found)
    {
      PetscStrToArray(data,'=',&j, &data2);
      if (j==2) {
        PetscStrToArray(data2[1],';',&k, &data3);
        PetscStrcpy(str_value_space, data3[0]); 
        PetscStrToArrayDestroy(k, data3);}
      else { PetscErrorPrintf("Error in input file. Size of variable '%s' is %d should be 2.\n", str_search, j);
             SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
      PetscStrToArrayDestroy(j, data2);
      break;
    }
  }
  ierr = PetscFClose(PETSC_COMM_WORLD, file); CHKERRQ(ierr);

  if PetscNot(found){ PetscErrorPrintf("Error in input file. Variable '%s' not found.\n", str_search);
                      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }

  PetscStrrchr(str_value_space,' ',&str_value_temp); // Removes first space
  if (BoolLower)  PetscStrtolower(str_value_temp);

  PetscStrcpy(str_value, str_value_temp);
  return 0;
}


PetscErrorCode SearchInputAsciiFileIfExists(const char* file_name, const char str_search[], char str_value[],const PetscBool BoolLower)
{
  /* Searches for string in a file. If BoolLower is true ==> change to lower case. */
  PetscErrorCode ierr = 0;
  PetscInt        i;
  PetscBool   found = PETSC_FALSE;
  char data[256]="";
  char **data2;
  char **data3;
  PetscInt j, k;
  char str_value_space[256]="";
  char *str_value_temp;

  FILE *file;
  ierr = PetscFOpen(PETSC_COMM_WORLD, file_name,"r", &file); CHKERRQ(ierr);
  for (i=0; i < 1000; i++)
  {
    PetscSynchronizedFGets(PETSC_COMM_WORLD, file, 256, data);
    PetscStrInList(str_search, data,' ', &found);
    if (found)
    {
      PetscStrToArray(data,'=',&j, &data2);
      if (j==2) {
        PetscStrToArray(data2[1],';',&k, &data3);
        PetscStrcpy(str_value_space, data3[0]); 
        PetscStrToArrayDestroy(k, data3);}
      else { PetscErrorPrintf("Error in input file. Size of variable '%s' is %d should be 2.\n", str_search, j);
             SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
      PetscStrToArrayDestroy(j, data2);
      break;
    }
  }
  ierr = PetscFClose(PETSC_COMM_WORLD, file); CHKERRQ(ierr);

  if (found)
  {
    PetscStrrchr(str_value_space,' ',&str_value_temp); // Removes first space
    if (BoolLower)  PetscStrtolower(str_value_temp);
    PetscStrcpy(str_value, str_value_temp);
    return 0;
  }
  else
    return 1;
}

PetscInt SearchInputAsciiFileForNumberOfOccurances(const char* file_name, const char str_search[])
{
  /* Searches for string in a file. If BoolLower is true ==> change to lower case. */
  PetscErrorCode ierr = 0;
  PetscInt        i;
  PetscBool   found = PETSC_FALSE;
  char data[256]="";
  char **data2;
  char **data3;
  PetscInt j, k;
  PetscInt N = 0;

  FILE *file;

  ierr = PetscFOpen(PETSC_COMM_WORLD, file_name,"r", &file); CHKERRQ(ierr);
  for (i=0; i < 1000; i++)
  {
    PetscSynchronizedFGets(PETSC_COMM_WORLD, file, 256, data);
    PetscStrInList(str_search, data,' ', &found);
    if (found)
    {
      PetscStrToArray(data,'=',&j, &data2);
      if (j==2) 
      {
        PetscStrToArray(data2[1],';',&k, &data3);
        PetscStrToArrayDestroy(k, data3);
      }
      else { PetscErrorPrintf("Error in input file. Size of variable '%s' with '=' is %d should be 2.\n", str_search, j);
             SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
      PetscStrToArrayDestroy(j, data2);
      N++;
      found=PETSC_FALSE;
    }
  }
  ierr = PetscFClose(PETSC_COMM_WORLD, file); CHKERRQ(ierr);

  if (N==0){ PetscErrorPrintf("Error in input file. Variable '%s' not found.\n", str_search);
             SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }

  return N;
}

PetscErrorCode SearchInputAsciiFileForNthOccurance(const char* file_name, const char str_search[], char str_value[], const PetscInt n)
{  
  PetscErrorCode  ierr = 0;
  PetscInt        i;
  PetscBool   found = PETSC_FALSE;
  char data[256]="";
  char **data2;
  char **data3;
  PetscInt j, k;
  char str_value_space[256]="";
  char *str_value_temp;
  PetscInt N = 0;

  FILE *file;
  ierr = PetscFOpen(PETSC_COMM_WORLD, file_name,"r", &file); CHKERRQ(ierr);
  for (i=0; i < 1000; i++)
  {
    PetscSynchronizedFGets(PETSC_COMM_WORLD, file, 256, data);
    PetscStrInList(str_search, data,' ', &found);
    if (found)
    {
      PetscStrToArray(data,'=',&j, &data2);
      if (n==N)
      {
        if (j==2) 
        {
          PetscStrToArray(data2[1],';',&k, &data3);
          PetscStrcpy(str_value_space, data3[0]); 
          PetscStrToArrayDestroy(k, data3);
          PetscStrToArrayDestroy(j, data2);
          break;
        }
        else { PetscErrorPrintf("Error in input file. Size of variable '%s' with '=' is %d should be 2.\n", str_search, j);
               SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }}
      PetscStrToArrayDestroy(j, data2);
      N++;
      found=PETSC_FALSE;
    }
  }
  ierr = PetscFClose(PETSC_COMM_WORLD, file); CHKERRQ(ierr);

  if (n!=N){ PetscErrorPrintf("Error in input file. The %d occrurance of variable '%s' not found.\n", N, str_search);
                      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }

  PetscStrrchr(str_value_space,' ',&str_value_temp); // Removes first space
  PetscStrtolower(str_value_temp);

  PetscStrcpy(str_value, str_value_temp);
  return 0;
}


PetscErrorCode SearchInputAsciiFileForStr(const char* file_name, const char str_search[], char str_value[],const PetscBool BoolLower)
{
  /* Searches for string in a file. If BoolLower is true ==> change to lower case. */
  SearchInputAsciiFile(file_name, str_search, str_value, BoolLower);
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, str_value); 
  return 0;
}

PetscErrorCode SearchInputAsciiFileForInt(const char* file_name, const char str_search[], PetscInt *value)
{
  /* Searches for a PetscInt in a file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  SearchInputAsciiFile(file_name, str_search, str_value, PETSC_FALSE);
  ierr = PetscOptionsStringToInt(str_value, value); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %d\n", str_search, *value); 
  return 0;
}


PetscErrorCode SearchInputAsciiFileForIntIfFound(const char* file_name, const char str_search[], PetscInt *value)
{
  /* Searches for a PetscBool in a file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_FALSE);
  if (ierr == 0) { // If str found
    ierr = PetscOptionsStringToInt(str_value, value); CHKERRQ(ierr);
    PetscPrintf(PETSC_COMM_WORLD,"    %s: %d\n", str_search, *value); 
  }
  else {  // If str found Not found
    *value = -1;
  }
  return 0;
}

PetscErrorCode SearchInputAsciiFileForScalar(const char* file_name, const char str_search[], PetscScalar *value)
{
  /* Searches for a PetscScalar in a file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  SearchInputAsciiFile(file_name, str_search, str_value, PETSC_FALSE);
  ierr = PetscOptionsStringToScalar(str_value, value); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %lf\n", str_search, *value); 
  return 0;
}

PetscErrorCode SearchInputAsciiFileForBool(const char* file_name, const char str_search[], PetscBool *value)
{
  /* Searches for a PetscBool in a file.*/
  char str_value[256]="";
  PetscBool  flg;
  SearchInputAsciiFile(file_name, str_search, str_value, PETSC_TRUE);
  PetscStrcmp("true",str_value,&flg);
  if (flg) *value = PETSC_TRUE;
  else{
    PetscStrcmp("false",str_value,&flg);
    if (flg) *value = PETSC_FALSE;
    else {
      SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT , "Error! '%s' boolean type could not be determined. Use 'false' or 'true'", str_search);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT , "Error in input file");}
  }
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %d\n", str_search, *value); 
  return 0;
}

PetscErrorCode SearchInputAsciiFileForBoolIfFound(const char* file_name, const char str_search[], PetscBool *value)
{
  /* Searches for a PetscBool in a file.*/
  char str_value[256]="";
  PetscBool  flg;
  PetscErrorCode ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_TRUE);
  if (ierr == 0) { // If str found
    PetscStrcmp("true",str_value,&flg);
    if (flg) *value = PETSC_TRUE;
    else{
      PetscStrcmp("false",str_value,&flg);
      if (flg) *value = PETSC_FALSE;
      else {
        SETERRQ1(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT , "Error! %s boolean type could not be determined. Use 'false' or 'true'", str_search);
        SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT , "Error in input file");}
    }
  }
  else {  // If str found Not found
    *value = PETSC_FALSE;
  }
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, *value ? "true" : "false"); 
  return 0;
}



PetscErrorCode SearchInputAsciiFileFor3Scalar(const char* file_name, const char str_search[], PetscScalar scalar_arr[])
{
  /* Searches for an Enum MESH_TYPE in a file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  char **data;
  PetscInt j;
  SearchInputAsciiFile(file_name, str_search, str_value, PETSC_FALSE);
  PetscStrToArray(str_value,',',&j, &data);
  if (j==3) {
    ierr = PetscOptionsStringToScalar(data[0], &scalar_arr[0]); CHKERRQ(ierr);
    ierr = PetscOptionsStringToScalar(data[1], &scalar_arr[1]); CHKERRQ(ierr);
    ierr = PetscOptionsStringToScalar(data[2], &scalar_arr[2]); CHKERRQ(ierr);
  }
  else { PetscErrorPrintf("Error in input file. Size of variable '%s' is %d should be 3.\n", str_search, j);
         SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
  PetscStrToArrayDestroy(j, data);
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %lf, %lf, %lf\n", str_search, scalar_arr[0], scalar_arr[1], scalar_arr[2]); 
  return 0;
}

PetscErrorCode SearchInputAsciiFileForScalarDataFile(const char* file_name, const char str_search[], System_Struct *system, PetscInt *Data_size, PetscScalar **Data)
{
  /* Searches for a PetscInt in a file.*/
  char data_file_name[256]="";
  PetscScalar neg_one = -1.0;
  SearchInputAsciiFile(file_name, str_search, data_file_name, PETSC_FALSE);
  
  *Data_size = FindSizeOfDataFile(data_file_name, system);
  PetscMalloc1(*Data_size, *&Data);
  for (PetscInt i=0; i<*Data_size; i++) *(Data[0]+i) = neg_one;
  ReadScalarDataFromFile(data_file_name, *Data_size, *Data, system);

  PetscPrintf(PETSC_COMM_WORLD,"    %s file read with size %d\n", data_file_name, *Data_size); 
  return 0;
}
PetscErrorCode SearchInputAsciiFileForEnumFuncTypeIfFound(const char* file_name, const char str_search[], FUNCTION_TYPE *value)
{
  /* Searches for an Enum VC Force type in file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  PetscBool flg_c1, flg_c2, flg_l1, flg_e1, flg_e2; 

  ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_TRUE);
  if (ierr == 0) // If str found
  {
    PetscStrcmp(str_value, "const",      &flg_c1);
    PetscStrcmp(str_value, "constant",   &flg_c2);

    PetscStrcmp(str_value, "linear",      &flg_l1);

    PetscStrcmp(str_value, "exp",         &flg_e1);
    PetscStrcmp(str_value, "exponential", &flg_e2);

    if      (flg_c1 || flg_c2)    *value = FUNC_CONSTANT;    
    else if (flg_l1)              *value = FUNC_LINEAR;
    else if (flg_e1 || flg_e2)    *value = FUNC_EXP;
    else{ 
      PetscErrorPrintf("Error in input file. Incorrect Enum FUNCTION_TYPE '%s' for variable '%s'.\n The correct inputs are: 'constant', 'linear',  or 'exponential'", str_value, str_search);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
  }
  else 
  {  // If str found Not found
    *value = FUNC_CONSTANT;
  }
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, FUNCTION_TYPE_CHAR[*value]); 

  return 0;
}



PetscErrorCode SearchInputAsciiFileForEnumSystemType(const char* file_name, const char str_search[], SYSTEM_TYPE *value)
{
  /* Searches for an Enum SYSTEM_TYPE in a file.*/
  char str_value[256]="";
  PetscBool flg_ibm, flg_stokes, flg_diff, flg_le;
  SearchInputAsciiFile(file_name, str_search, str_value, PETSC_TRUE);

  PetscStrcmp(str_value, "ibm",    &flg_ibm);
  PetscStrcmp(str_value, "stokes", &flg_stokes);
  PetscStrcmp(str_value, "diff",   &flg_diff);
  PetscStrcmp(str_value, "linear_elasticity",   &flg_le);
  
  if      (flg_ibm)     *value = IBM;    
  else if (flg_stokes)  *value = STOKES;
  else if (flg_diff)    *value = DIFF;
  else if (flg_le)      *value = LINEAR_ELASTICITY;
  else { 
    PetscErrorPrintf("Error in input file. Incorrect Enum SYSTEM_TYPE '%s' for variable '%s'.\n", str_value, str_search);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }

  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, SYSTEM_TYPE_CHAR[*value]); 
  return 0;
}

PetscErrorCode SearchInputAsciiFileForEnumMeshType(const char* file_name, const char str_search[], MESH_TYPE *value)
{
  /* Searches for an Enum MESH_TYPE in a file.*/
  char str_value[256]="";
  PetscBool flg_gmsh, flg_vtk, flg_foam, flg_gmsh_old;
  SearchInputAsciiFile(file_name, str_search, str_value, PETSC_TRUE);
  
  PetscStrcmp(str_value, "gmsh", &flg_gmsh);
  PetscStrcmp(str_value, "vtk",  &flg_vtk);
  PetscStrcmp(str_value, "foam", &flg_foam);
  PetscStrcmp(str_value, "gmsh_old", &flg_gmsh_old);

  if      (flg_gmsh)      *value = GMSH;    
  else if (flg_vtk)       *value = VTK;
  else if (flg_foam)      *value = FOAM;
  else if (flg_gmsh_old)  *value = GMSH_OLD;
  else{ 
    PetscErrorPrintf("Error in input file. Incorrect Enum MESH_TYPE '%s' for variable '%s'.\n", str_value, str_search);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }

  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, MESH_TYPE_CHAR[*value]); 
  return 0;
}



// Fluid Enumerate Functions
//=================================================================
PetscErrorCode SearchInputAsciiFileForEnumFLSolverType(const char* file_name, const char str_search[], FLUID_SOLVER_TYPE *value)
{
  /* Searches for an Enum FLUID_SOLVER_TYPE in a file.*/
  char str_value[256]="";
  PetscBool flg_imp, flg_exp_lump, flg_imp_exp_lump, flg_exp;
  SearchInputAsciiFile(file_name, str_search, str_value, PETSC_TRUE);

  PetscStrcmp(str_value, "implicit",                  &flg_imp);
  PetscStrcmp(str_value, "explicit_lumped",           &flg_exp_lump);
  PetscStrcmp(str_value, "implicit_explicit_lumped",  &flg_imp_exp_lump);
  PetscStrcmp(str_value, "explicit",                  &flg_exp);
  
  if      (flg_imp)           *value = IMPLICIT;    
  else if (flg_exp_lump)      *value = EXPLICIT_LUMPED;
  else if (flg_imp_exp_lump)  *value = IMPLICIT_EXPLICIT_LUMPED;
  else if (flg_exp)           *value = EXPLICIT;
  else{ 
    PetscErrorPrintf("Error in input file. Incorrect Enum FLUID_SOLVER_TYPE '%s' for variable '%s'.\n", str_value, str_search);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }

  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, FLUID_SOLVER_TYPE_CHAR[*value]); 
  return 0;
}


PetscErrorCode SearchInputAsciiFileForEnumFLSolverMethod(const char* file_name, const char str_search[], FLUID_SOLVER_METHOD *value)
{
  /* Searches for an Enum FLUID_SOLVER_METHOD in a file.*/
  char str_value[256]="";
  PetscBool flg_direct, flg_nest, flg_mono, flg_mono_pre;
  SearchInputAsciiFile(file_name, str_search, str_value, PETSC_TRUE);
  
  PetscStrcmp(str_value, "direct", &flg_direct);
  PetscStrcmp(str_value, "nest",   &flg_nest);
  PetscStrcmp(str_value, "mono",   &flg_mono);
  PetscStrcmp(str_value, "mono_pre",   &flg_mono_pre);

  if      (flg_direct)   *value = DIRECT;    
  else if (flg_nest)     *value = NEST;
  else if (flg_mono)     *value = MONO;
  else if (flg_mono_pre) *value = MONO_PRE;
  else{  
    PetscErrorPrintf("Error in input file. Incorrect Enum FLUID_SOLVER_METHOD '%s' for variable '%s'.\n", str_value, str_search);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }

  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, FLUID_SOLVER_METHOD_CHAR[*value]); 
  return 0;
}


PetscErrorCode SearchInputAsciiFileForEnumFLElemType(const char* file_name, const char str_search[], FLUID_ELEMENT_TYPE *value)
{
  /* Searches for an Enum MESH_TYPE in a file.*/
  char str_value[256]="";
  PetscBool flg_p0, flg_p1, flg_p1b, flg_p2;
  SearchInputAsciiFile(file_name, str_search, str_value, PETSC_TRUE);
  
  PetscStrcmp(str_value, "p0",  &flg_p0);
  PetscStrcmp(str_value, "p1",  &flg_p1);
  PetscStrcmp(str_value, "p1b", &flg_p1b);
  PetscStrcmp(str_value, "p2",  &flg_p2);
  
  if      (flg_p0)  *value = P0;    
  else if (flg_p1)  *value = P1;
  else if (flg_p1b) *value = P1b;
  else if (flg_p2)  *value = P2;
  else{ 
    PetscErrorPrintf("Error in input file. Incorrect Enum ELEMENT_TYPE '%s' for variable '%s'.\n", str_value, str_search);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }

  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, FLUID_ELEMENT_TYPE_CHAR[*value]); 
  return 0;
}


PetscErrorCode SearchInputAsciiFileForEnumFLViscMethod(const char* file_name, const char str_search[], FLUID_VISCOSITY_METHOD *value)
{
  /* Searches for an Enum FLUID_SOLVER_METHOD in a file.*/
  char str_value[256]="";
  PetscBool flg_1c, flg_2c1e, flg_2c2e; 

  SearchInputAsciiFile(file_name, str_search, str_value, PETSC_TRUE);

  PetscStrcmp(str_value, "one_constant",                    &flg_1c);
  PetscStrcmp(str_value, "two_constants_one_ellipsoid",     &flg_2c1e);
  PetscStrcmp(str_value, "two_constants_two_ellipsoids",    &flg_2c2e);

  if      (flg_1c)    *value = ONE_CONSTANT;    
  else if (flg_2c1e)  *value = TWO_CONSTANTS_ONE_ELLIPSOID;
  else if (flg_2c2e)  *value = TWO_CONSTANTS_TWO_ELLIPSOIDS;
  else{ 
    PetscErrorPrintf("Error in input file. Incorrect Enum FLUID_VISCOSITY_METHOD_CHAR '%s' for variable '%s'.\n The correct inputs are: 'one_constant', 'two_constants_one_ellipsoid', 'two_constants_two_ellipsoids'", str_value, str_search);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");   }

  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, FLUID_VISCOSITY_METHOD_CHAR[*value]); 
  return 0;
}


// Solid Enumerate Functions
//=================================================================
PetscErrorCode SearchInputAsciiFileForEnumSDForceLETypeIfFound(const char* file_name, const char str_search[], SOLID_FORCE_LE_TYPE *value)
{
  /* Searches for an Enum LE Force type in file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  PetscBool flg_01, flg_02, flg_03, flg_04, flg_le_fd, flg_le_id, flg_le_imk, flg_le_id_cond1, flg_le_id_cond2; 
  PetscBool flg_le_fd2, flg_le_id2, flg_le_imk2, flg_le_id_cond3, flg_le_id_cond4; 
  PetscBool flg_le_fd3, flg_le_id3, flg_le_imk3; 

  ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_TRUE);
  if (ierr == 0) // If str found
  {
    PetscStrcmp(str_value, "f0",      &flg_01);
    PetscStrcmp(str_value, "null",    &flg_02);
    PetscStrcmp(str_value, "false",   &flg_03);
    PetscStrcmp(str_value, "no",      &flg_04);

    PetscStrcmp(str_value, "fixed_disp",   &flg_le_fd);;
    PetscStrcmp(str_value, "le_fixed_disp",   &flg_le_fd2);
    PetscStrcmp(str_value, "le_fd",   &flg_le_fd3);

    PetscStrcmp(str_value, "initial_disp",   &flg_le_id);
    PetscStrcmp(str_value, "le_initial_disp",   &flg_le_id2);
    PetscStrcmp(str_value, "le_id",   &flg_le_id3);

    PetscStrcmp(str_value, "initial_disp_multi_k", &flg_le_imk);
    PetscStrcmp(str_value, "le_initial_disp_multi_k", &flg_le_imk2);
    PetscStrcmp(str_value, "le_imk", &flg_le_imk3);

    PetscStrcmp(str_value, "initial_disp_conditional", &flg_le_id_cond1);
    PetscStrcmp(str_value, "le_initial_disp_conditional", &flg_le_id_cond3);
    PetscStrcmp(str_value, "initial_disp_cond", &flg_le_id_cond2);
    PetscStrcmp(str_value, "le_initial_disp_cond", &flg_le_id_cond4);

    if      (flg_01 || flg_02 || flg_03 || flg_04)      *value = FORCE_LE_NULL;    
    else if (flg_le_fd || flg_le_fd2 || flg_le_fd3)     *value = FORCE_LE_FIXED_DISP;
    else if (flg_le_id || flg_le_id2 || flg_le_id3)     *value = FORCE_LE_INITIAL_DISP;
    else if (flg_le_imk || flg_le_imk2 || flg_le_imk3)  *value = FORCE_LE_INITIAL_DISP_MULTI;
    else if (flg_le_id_cond1 || flg_le_id_cond2 || flg_le_id_cond3 || flg_le_id_cond4) *value = FORCE_LE_INITIAL_DISP_COND;
    else{ 
      PetscErrorPrintf("Error in input file. Incorrect Enum SOLID_FORCE_LE_TYPE '%s' for variable '%s'.\n The correct inputs are: 'false', 'fixed_disp', 'initial_disp', 'initial_disp_multi_k', and 'initial_disp_conditional'", str_value, str_search);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
  }
  else 
  {  // If str found Not found
    *value = FORCE_LE_NULL;
  }
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, SOLID_FORCE_LE_CHAR[*value]); 
  return 0;
}

PetscErrorCode SearchInputAsciiFileForEnumSDForceSTTypeIfFound(const char* file_name, const char str_search[], SOLID_FORCE_ST_TYPE *value)
{
  /* Searches for an Enum ST Force type in file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  PetscBool flg_st_no1, flg_st_no2, flg_st_c1, flg_st_c2, flg_st_c3, flg_st_m1, flg_st_m2, flg_st_m3;
  PetscBool flg_st_cp1, flg_st_cp2, flg_st_cp3, flg_st_mp1, flg_st_mp2, flg_st_mp3; 

  ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_TRUE);
  if (ierr == 0) // If str found
  {
    PetscStrcmp(str_value, "false",      &flg_st_no1);
    PetscStrcmp(str_value, "no_st",      &flg_st_no2);

    PetscStrcmp(str_value, "st_const",   &flg_st_c1);
    PetscStrcmp(str_value, "const",      &flg_st_c2);
    PetscStrcmp(str_value, "constant",   &flg_st_c3);

    PetscStrcmp(str_value, "st_multi",   &flg_st_m1);
    PetscStrcmp(str_value, "multi",      &flg_st_m2);
    PetscStrcmp(str_value, "multiple",   &flg_st_m3);

    PetscStrcmp(str_value, "st_const_perc",         &flg_st_cp1);
    PetscStrcmp(str_value, "const_perc",            &flg_st_cp2);
    PetscStrcmp(str_value, "constant_percentage",   &flg_st_cp3);

    PetscStrcmp(str_value, "st_multi_perc",   &flg_st_mp1);
    PetscStrcmp(str_value, "multi_perc",      &flg_st_mp2);
    PetscStrcmp(str_value, "multiple_percentage",   &flg_st_mp3);

    if      (flg_st_no1 || flg_st_no2)                  *value = FORCE_ST_NULL;    
    else if (flg_st_c1  || flg_st_c2  || flg_st_c3)     *value = FORCE_ST_CONST;
    else if (flg_st_m1  || flg_st_m2  || flg_st_m3)     *value = FORCE_ST_MULTI; 
    else if (flg_st_cp1 || flg_st_cp2 || flg_st_cp3)    *value = FORCE_ST_CONST_PERC;
    else if (flg_st_mp1 || flg_st_mp2 || flg_st_mp3)    *value = FORCE_ST_MULTI_PERC;
    else{ 
      PetscErrorPrintf("Error in input file. Incorrect Enum SOLID_FORCE_ST_TYPE '%s' for variable '%s'.\n The correct inputs are: 'false', 'no_st', 'st_c', 'st_m'", str_value, str_search);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
  }
  else 
  {  // If str found Not found
    *value = FORCE_ST_NULL;
  }
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, SOLID_FORCE_ST_CHAR[*value]); 

  return 0;
}



PetscErrorCode SearchInputAsciiFileForEnumSDForceVCTypeIfFound(const char* file_name, const char str_search[], SOLID_FORCE_VC_TYPE *value)
{
  /* Searches for an Enum VC Force type in file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  PetscBool flg_vc_no1, flg_vc_no2, flg_vc_c1, flg_vc_c2, flg_vc_c3, flg_vc_c4; 

  ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_TRUE);
  if (ierr == 0) // If str found
  {
    PetscStrcmp(str_value, "false",      &flg_vc_no1);
    PetscStrcmp(str_value, "no_vc",      &flg_vc_no2);

    PetscStrcmp(str_value, "vc_const",   &flg_vc_c1);
    PetscStrcmp(str_value, "const",      &flg_vc_c2);
    PetscStrcmp(str_value, "constant",   &flg_vc_c3);
    PetscStrcmp(str_value, "true",       &flg_vc_c4);

    if      (flg_vc_no1 || flg_vc_no2)                           *value = FORCE_VC_NULL;    
    else if (flg_vc_c1 || flg_vc_c2 || flg_vc_c3 || flg_vc_c4)   *value = FORCE_VC_CONST;
    else{ 
      PetscErrorPrintf("Error in input file. Incorrect Enum SOLID_FORCE_VC_TYPE '%s' for variable '%s'.\n The correct inputs are: 'false', 'no_vc', 'st_c', 'st_m'", str_value, str_search);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
  }
  else 
  {  // If str found Not found
    *value = FORCE_VC_NULL;
  }
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, SOLID_FORCE_VC_CHAR[*value]); 

  return 0;
}


PetscErrorCode SearchInputAsciiFileForEnumSDForceBFTypeIfFound(const char* file_name, const char str_search[], SOLID_FORCE_BF_TYPE *value)
{
  /* Searches for an Enum VC Force type in file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  PetscBool flg_bf_no1, flg_bf_no2, flg_bf_a1, flg_bf_a2, flg_bf_bp1, flg_bf_bp2; 

  ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_TRUE);
  if (ierr == 0) // If str found
  {
    PetscStrcmp(str_value, "false",    &flg_bf_no1);
    PetscStrcmp(str_value, "no_bf",    &flg_bf_no2);

    PetscStrcmp(str_value, "all", &flg_bf_a1);
    PetscStrcmp(str_value, "all_points",    &flg_bf_a2);

    PetscStrcmp(str_value, "boundary_points",    &flg_bf_bp1);
    PetscStrcmp(str_value, "boundary_pts", &flg_bf_bp2);

    if      (flg_bf_no1 || flg_bf_no2)    *value = FORCE_BF_NULL;    
    else if (flg_bf_a1 || flg_bf_a2)      *value = FORCE_BF_ALL_POINTS;
    else if (flg_bf_bp1 || flg_bf_bp2)    *value = FORCE_BF_BOUNDARY_POINTS;
    else{ 
      PetscErrorPrintf("Error in input file. Incorrect Enum SOLID_FORCE_BF_TYPE '%s' for variable '%s'.\n The correct inputs are: 'false', 'const', or 'boundary_pts", str_value, str_search);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
  }
  else 
  {  // If str found Not found
    *value = FORCE_BF_NULL;
  }
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, SOLID_FORCE_BF_CHAR[*value]); 

  return 0;
}



PetscErrorCode SearchInputAsciiFileForEnumSDForceSFTypeIfFound(const char* file_name, const char str_search[], SOLID_FORCE_SF_TYPE *value)
{
  /* Searches for an Enum SF Force type in file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  PetscBool flg_sf_no1, flg_sf_no2, flg_sf_e, flg_sf_2e; 

  ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_TRUE);
  if (ierr == 0) // If str found
  {
    PetscStrcmp(str_value, "false",      &flg_sf_no1);
    PetscStrcmp(str_value, "no_sf",      &flg_sf_no2);

    PetscStrcmp(str_value, "one_ellipsoid",   &flg_sf_e);
    PetscStrcmp(str_value, "two_ellipsoids",  &flg_sf_2e);

    if      (flg_sf_no1 || flg_sf_no2)    *value = FORCE_SF_NULL;    
    else if (flg_sf_e)                    *value = FORCE_SF_ELLIPSOID;
    else if (flg_sf_2e)                   *value = FORCE_SF_TWO_ELLIPSOIDS;
    else{ 
      PetscErrorPrintf("Error in input file. Incorrect Enum SOLID_FORCE_SF_TYPE '%s' for variable '%s'.\n The correct inputs are: 'false', 'no_sf', 'st_c', 'st_m'", str_value, str_search);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
  }
  else 
  {  // If str found Not found
    *value = FORCE_SF_NULL;
  }
  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, SOLID_FORCE_SF_CHAR[*value]); 

  return 0;
}



PetscErrorCode SearchInputAsciiFileForEnumSDConstrainTypeIfFound(const char* file_name, const char str_search[], SOLID_CONSTRAIN_TYPE *value)
{
  /* Searches for an Enum FLUID_SOLVER_METHOD in a file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  PetscBool flg_false, flg_nc, flg_1e, flg_2e, flg_axs; 

  ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_TRUE);
  if (ierr == 0) // If str found
  {
    PetscStrcmp(str_value, "false",              &flg_false);
    PetscStrcmp(str_value, "no_constrain",       &flg_nc);
    PetscStrcmp(str_value, "one_ellipsoid",      &flg_1e);
    PetscStrcmp(str_value, "two_ellipsoids",     &flg_2e);
    PetscStrcmp(str_value, "axis_constrain",     &flg_axs);

    if      (flg_nc || flg_false)     *value = CONSTRAIN_NULL;    
    else if (flg_1e)                  *value = CONSTRAIN_ELLIPSOID;
    else if (flg_2e)                  *value = CONSTRAIN_TWO_ELLIPSOIDS;
    else if (flg_axs)                 *value = CONSTRAIN_AXIS;
    else{ 
      PetscErrorPrintf("Error in input file. Incorrect Enum SOLID_CONSTRAIN_TYPE '%s' for variable '%s'.\n The correct inputs are: 'no_constrain', 'one_ellipsoid', 'two_ellipsoids' and 'axis_constrain'", str_value, str_search);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
  }
  else 
  {  // If str found Not found
    *value = CONSTRAIN_NULL;
  }

  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, SOLID_CONSTRAIN_CHAR[*value]); 
  return 0;
}

PetscErrorCode SearchInputAsciiFileForEnumSDLNodeTypeIfFound(const char* file_name, const char str_search[], SOLID_LNODE_DYNAMICS_TYPE *value)
{
  /* Searches for an Enum FLUID_SOLVER_METHOD in a file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  PetscBool flg_false, flg_n, flg_a, flg_c1, flg_c2; 

  ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_TRUE);
  if (ierr == 0) // If str found
  {
    PetscStrcmp(str_value, "false",         &flg_false);
    PetscStrcmp(str_value, "no",            &flg_n);
    PetscStrcmp(str_value, "always",        &flg_a);
    PetscStrcmp(str_value, "conditional",   &flg_c1);
    PetscStrcmp(str_value, "cond",          &flg_c2);

    if      (flg_n || flg_false)   *value = LNODE_DYNAMICS_NULL;    
    else if (flg_a)                *value = LNODE_DYNAMICS_ALWAYS;
    else if (flg_c1 || flg_c2)     *value = LNODE_DYNAMICS_COND;
    else{ 
      PetscErrorPrintf("Error in input file. Incorrect Enum SOLID_LNODE_DYNAMICS_TYPE '%s' for variable '%s'.\n The correct inputs are: 'false', 'always', and 'conditional'", str_value, str_search);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
  }
  else 
  {  // If str found Not found
    *value = LNODE_DYNAMICS_NULL;
  }

  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, SOLID_LNODE_DYNAMICS_CHAR[*value]); 
  return 0;
}

PetscErrorCode SearchInputAsciiFileForEnumSDNoiseTypeIfFound(const char* file_name, const char str_search[], SOLID_NOISE_TYPE *value)
{
  /* Searches for an Enum FLUID_SOLVER_METHOD in a file.*/
  PetscErrorCode ierr = 0;  
  char str_value[256]="";
  PetscBool flg_false, flg_n, flg_c0, flg_c1, flg_c2; 

  ierr = SearchInputAsciiFileIfExists(file_name, str_search, str_value, PETSC_TRUE);
  if (ierr == 0) // If str found
  {
    PetscStrcmp(str_value, "false",         &flg_false);
    PetscStrcmp(str_value, "no",            &flg_n);
    PetscStrcmp(str_value, "true",          &flg_c0);
    PetscStrcmp(str_value, "conditional",   &flg_c1);
    PetscStrcmp(str_value, "cond",          &flg_c2);

    if      (flg_n || flg_false)          *value = NOISE_NULL;    
    else if (flg_c0 || flg_c1 || flg_c2)  *value = NOISE_CONST;
    else{ 
      PetscErrorPrintf("Error in input file. Incorrect Enum SOLID_NOISE_TYPE '%s' for variable '%s'.\n The correct inputs are: 'false', or 'true'", str_value, str_search);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_USER_INPUT,"Error in input file.");  }
  }
  else 
  {  // If str found Not found
    *value = NOISE_NULL;
  }

  PetscPrintf(PETSC_COMM_WORLD,"    %s: %s\n", str_search, SOLID_NOISE_CHAR[*value]); 
  return 0;
}

#endif
