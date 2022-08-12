#ifndef VTK_WriteSolidFiles
#define VTK_WriteSolidFiles


PetscErrorCode VTKGetFluidElemInfo(PetscInt *VTKInternalElements, PetscInt *VTKBoundaryElements, 
    PetscInt *InternalSize, PetscInt *BoundarySize, const PetscInt Dim, const FLUID_ELEMENT_TYPE ET)
{
  
  if (Dim == 2 && ET == P1) 
    *VTKInternalElements = 5,  *VTKBoundaryElements=3, *InternalSize=3, *BoundarySize=2;
  else if (Dim == 2 && ET == Q1) 
    *VTKInternalElements = 7,  *VTKBoundaryElements=4, *InternalSize=4, *BoundarySize=2;
  else if (Dim == 3 && ET == P1) 
    *VTKInternalElements = 10, *VTKBoundaryElements=5, *InternalSize=4, *BoundarySize=3;
  else if (Dim == 3 && ET == Q1) 
    *VTKInternalElements = 12, *VTKBoundaryElements=9, *InternalSize=6, *BoundarySize=2;
  else if (Dim == 2 && ET == P2)
    *VTKInternalElements = 22, *VTKBoundaryElements=21, *InternalSize=6, *BoundarySize=3;
  else if (Dim == 3 && ET == P2)
    *VTKInternalElements = 24, *VTKBoundaryElements=22, *InternalSize=10, *BoundarySize=6;
  else
    { PetscErrorPrintf("The VTK element (Dim = %d, Type = %s) is not defined.\n", Dim, FLUID_ELEMENT_TYPE_CHAR[ET]);
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT,"Error in GetVTKElemInfo function.\n"); }
  return 0;
}

PetscErrorCode VTKGetSolidElemInfo(PetscInt *VTKElements, const SOLID_ELEMENT_TYPE ET)
{
  
  if (ET == POINT)          *VTKElements = 1;
  else if (ET == LINE)      *VTKElements = 3;
  else if (ET == TRIANGLE)  *VTKElements = 5;
  else  SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG, "The only cells that are defined are:\n\t1- LINE\n\t2-TRIANGLE\n");
  return 0;
}

void VTKWriteIntro(char *VTKFile_Name, const char *VTKFile_Title)
{
  FILE *fd;
  fd = fopen(VTKFile_Name,"w");
  fprintf(fd, "# vtk DataFile Version 2.0\n");
  fprintf(fd, "%s, Created by Mohamad Ibrahim Cheikh\n", VTKFile_Title);
  fprintf(fd, "ASCII\n");
  fprintf(fd, "DATASET UNSTRUCTURED_GRID\n");
  fclose(fd);
}

void VTKWriteIntroParallel(char *VTKFile_Name, const char *VTKFile_Title)
{
  PetscViewer    viewer;
  PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
  PetscViewerSetType(viewer, PETSCVIEWERASCII);
  PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
  PetscViewerFileSetName(viewer, VTKFile_Name);
  PetscViewerASCIIPrintf(viewer, "# vtk DataFile Version 2.0\n");
  PetscViewerASCIIPrintf(viewer, "%s, Created by Mohamad Ibrahim Cheikh\n", VTKFile_Title);
  PetscViewerASCIIPrintf(viewer, "ASCII\n");
  PetscViewerASCIIPrintf(viewer, "DATASET UNSTRUCTURED_GRID\n");
  PetscViewerDestroy(&viewer);
}


void VTKWriteNodes(char *VTKFile_Name,const PetscInt N, PetscScalar *x, PetscScalar *y, PetscScalar *z)
{
  FILE *fd;
  fd = fopen(VTKFile_Name,"a");
  fprintf(fd, "POINTS %d double\n", N);
  for(PetscInt i=0; i<N; i++)
    fprintf(fd, "%g %g %g\n", (double)x[i], (double)y[i], (double)z[i]);
  fprintf(fd, "\n");
  fclose(fd);
}


void VTKWriteNodesP(char *VTKFile_Name, PetscInt N, PetscScalar *x, PetscScalar *y, PetscScalar *z, PetscInt *LI_PtoU)
{
  /* Only for P nodes */
  FILE *fd;
  fd = fopen(VTKFile_Name,"a");
  fprintf(fd, "POINTS %d double\n", N);
  for(PetscInt i=0; i< N; i++) {
    const PetscInt LI = LI_PtoU[i];
    fprintf(fd, "%g %g %g\n", (double)x[LI], (double)y[LI], (double)z[LI]); }
  fprintf(fd, "\n");
  fclose(fd);
} 

void VTKWriteSolidElements(char *VTKFile_Name, Solid_Struct *solid)
{
  PetscInt N_Cells = solid->N_total_lns + solid->N_total_trs;
  FILE *fd;
  fd = fopen(VTKFile_Name,"a");

  fprintf(fd, "CELLS %d %d\n", N_Cells, solid->N_total_lns*3 + solid->N_total_trs*4);
  for(PetscInt l=0; l < solid->N_total_lns; l++)
    fprintf(fd, "2 %d %d\n",  solid->lns[l][0], solid->lns[l][1]);
  for(PetscInt t=0; t < solid->N_total_trs; t++)
    fprintf(fd, "3 %d %d %d\n",  solid->trs[t][0], solid->trs[t][1], solid->trs[t][2]);
  fprintf(fd, "\n");

  fprintf(fd, "CELL_TYPES %d\n", N_Cells);
  for(PetscInt l=0; l < solid->N_total_lns; l++)  fprintf(fd, "3\n");
  for(PetscInt t=0; t < solid->N_total_trs; t++)  fprintf(fd, "5\n");
  fprintf(fd, "\n");

  fclose(fd);
}


void VTKWriteFluidElements(char *VTKFile_Name, const PetscInt Dim, Nodes_Struct *nodes, Elements_Struct *elements, 
  const PetscInt variable_number)
{  
  FILE *fd;
  PetscInt S_inter, S_bound;
  FLUID_ELEMENT_TYPE ET;

  if (variable_number == 1)      
    S_inter = elements->P_S_inter, S_bound = elements->P_S_bound, ET = elements->P_elem_type;
  else if (variable_number == 2) 
    S_inter = elements->U_S_inter, S_bound = elements->U_S_bound, ET = elements->U_elem_type;

  fd   = fopen(VTKFile_Name,"a");

  
  if (ET != P1b)
    fprintf(fd, "CELLS %d %d\n", (elements->N_inter_loc+elements->N_bound_loc), 
      ((S_inter+1)*elements->N_inter_loc + (S_bound+1)*elements->N_bound_loc) );
  else  
    fprintf(fd, "CELLS %d %d\n", (elements->N_inter_loc*(Dim+1)+elements->N_bound_loc), 
      ((S_inter)*elements->N_inter_loc*(Dim+1) + (S_bound+1)*elements->N_bound_loc) );
  
  if (variable_number == 2) 
    for(PetscInt i=0; i<elements->N_bound_loc; i++) {
      fprintf(fd, "%d", S_bound );
      for(PetscInt j=0; j<S_bound; j++)
        fprintf(fd, " %d", elements->bound[i][j] );
      fprintf(fd, "\n" );}
  else if (variable_number == 1) 
    for(PetscInt i=0; i<elements->N_bound_loc; i++) { 
      fprintf(fd, "%d", S_bound );
      for(PetscInt j=0; j<S_bound; j++)
        fprintf(fd, " %d", nodes->LI_UtoP[elements->bound[i][j]] );
      fprintf(fd, "\n" );}

  
  if (ET == P1b)
  {
    if (Dim == 2)
      for(PetscInt i=0; i<elements->N_inter_loc; i++) 
      {
        fprintf(fd, "%d %d %d %d\n", S_inter-1, elements->inter[i][0], elements->inter[i][1], elements->inter[i][3]); 
        fprintf(fd, "%d %d %d %d\n", S_inter-1, elements->inter[i][1], elements->inter[i][2], elements->inter[i][3]); 
        fprintf(fd, "%d %d %d %d\n", S_inter-1, elements->inter[i][2], elements->inter[i][0], elements->inter[i][3]); 
      }
    else if (Dim == 3)
      for(PetscInt i=0; i<elements->N_inter_loc; i++) 
      {
        fprintf(fd, "%d %d %d %d %d\n", S_inter-1, elements->inter[i][0], elements->inter[i][1], elements->inter[i][2], elements->inter[i][4]); 
        fprintf(fd, "%d %d %d %d %d\n", S_inter-1, elements->inter[i][1], elements->inter[i][2], elements->inter[i][3], elements->inter[i][4]); 
        fprintf(fd, "%d %d %d %d %d\n", S_inter-1, elements->inter[i][2], elements->inter[i][3], elements->inter[i][0], elements->inter[i][4]); 
        fprintf(fd, "%d %d %d %d %d\n", S_inter-1, elements->inter[i][3], elements->inter[i][0], elements->inter[i][1], elements->inter[i][4]); 
      }
  }
  else if (variable_number == 2) 
    for(PetscInt i=0; i<elements->N_inter_loc; i++) 
    {
      fprintf(fd, "%d", S_inter);
      for(PetscInt j=0; j<S_inter; j++)
        fprintf(fd, " %d", elements->inter[i][j]);
      fprintf(fd, "\n" );
    }
  else if (variable_number == 1) 
    for(PetscInt i=0; i<elements->N_inter_loc; i++) 
    {
      fprintf(fd, "%d", S_inter);
      for(PetscInt j=0; j<S_inter; j++)
        fprintf(fd, " %d", nodes->LI_UtoP[elements->inter[i][j]]);
      fprintf(fd, "\n" );
    }
  fprintf(fd, "\n" );
  

  fclose(fd);

  
  PetscInt VTK_inter, VTK_bound, Size_inter, Size_bound;
  if (ET!=P1b) VTKGetFluidElemInfo(&VTK_inter, &VTK_bound, &Size_inter, &Size_bound, Dim, ET);
  else         VTKGetFluidElemInfo(&VTK_inter, &VTK_bound, &Size_inter, &Size_bound, Dim, P1);

  PetscInt N_int, N_bou;
  if (ET!=P1b)  N_int = elements->N_inter_loc,          N_bou = elements->N_bound_loc;
  else          N_int = elements->N_inter_loc*(Dim+1),  N_bou = elements->N_bound_loc;

  fd   = fopen(VTKFile_Name,"a");
  fprintf(fd, "CELL_TYPES %d\n", N_int+N_bou); 
  for(PetscInt i=0; i<N_bou; i++)
    fprintf(fd, "%d\n", VTK_bound); 
  for(PetscInt i=0; i<N_int; i++)
    fprintf(fd, "%d\n", VTK_inter); 
  fprintf(fd, "\n" );
  fclose(fd);
}


void VTKWritePointData_Vector(char *VTKFile_Name, const char VarName[], PetscInt N, PetscScalar *x, PetscScalar *y, PetscScalar *z, PetscBool AddHeader)
{
  FILE *fd;
  fd = fopen(VTKFile_Name,"a");

  if (AddHeader) fprintf(fd, "POINT_DATA %d\n", N);

  fprintf(fd, "VECTORS %s double\n", VarName);
  for(PetscInt i=0; i < N; i++)
    fprintf(fd, "%g %g %g\n", (double)x[i], (double)y[i], (double)z[i]);
  fprintf(fd, "\n");
  fclose(fd);
}

void VTKWritePointData_Scalars(char *VTKFile_Name, const char VarName[], PetscInt N, PetscScalar *s, PetscBool AddHeader)
{
  FILE *fd;
  fd = fopen(VTKFile_Name,"a");

  if (AddHeader) fprintf(fd, "POINT_DATA %d\n", N);

  fprintf(fd, "SCALARS %s double\nLOOKUP_TABLE default\n", VarName);
  for(PetscInt i=0; i < N; i++)
    fprintf(fd, "%g\n", (double)s[i]);
  fprintf(fd, "\n");
  fclose(fd);
}

void VTKWritePointData_Ints(char *VTKFile_Name, const char VarName[], PetscInt N, PetscInt *L, PetscBool AddHeader)
{
  FILE *fd;
  fd = fopen(VTKFile_Name,"a");

  if (AddHeader) fprintf(fd, "POINT_DATA %d\n", N);

  fprintf(fd, "SCALARS %s int 1\nLOOKUP_TABLE default\n", VarName);
  for(PetscInt i=0; i < N; i++)
    fprintf(fd, "%d\n", L[i]);
  fprintf(fd, "\n");
  fclose(fd);
}

void VTKWritePointDataP_Ints(char *VTKFile_Name, const char VarName[], PetscInt N, PetscInt *L, PetscInt *LI_PtoU, PetscBool AddHeader)
{
  FILE *fd;
  fd = fopen(VTKFile_Name,"a");

  if (AddHeader) fprintf(fd, "POINT_DATA %d\n", N);

  fprintf(fd, "SCALARS %s int 1\nLOOKUP_TABLE default\n", VarName);
  for(PetscInt i=0; i < N; i++){
    const PetscInt LI = LI_PtoU[i];
    fprintf(fd, "%d\n", L[LI]);
  }
  fprintf(fd, "\n");
  fclose(fd);
}

void VTKWriteNodesParallel(char *VTKFile_Name, const PetscMPIInt size, const PetscMPIInt rank, Nodes_Struct *nodes, const PetscInt nodes_number)
{
  FILE *fd;
  PetscInt N_total_glo, N_total_loc;

  if(nodes_number == 1)       
    N_total_glo = nodes->N_P_total_glo, N_total_loc = nodes->N_P_total_loc;
  else if(nodes_number == 2)  
    N_total_glo = nodes->N_U_total_glo, N_total_loc = nodes->N_U_total_loc;

  for (PetscInt r=0; r< size; ++r) 
  {
    if (rank == r) 
    {
      fd = fopen(VTKFile_Name,"a");
      if (rank == 0) fprintf(fd, "POINTS %d double\n", N_total_glo);
      
      if (nodes_number == 2)  
      { 
        for(PetscInt i=0; i< N_total_loc; i++)
          if (nodes->t[i] < 2) 
            fprintf(fd, "%g %g %g\n", (double)nodes->x[i], (double)nodes->y[i], (double)nodes->z[i]);
      }
      else if (nodes_number == 1)  
      { 
        for(PetscInt i=0; i< N_total_loc; i++)
        { 
          const PetscInt i_P = nodes->LI_PtoU[i];
          if (nodes->t[i_P] < 2) 
            fprintf(fd, "%g %g %g\n", (double)nodes->x[i_P], (double)nodes->y[i_P], (double)nodes->z[i_P]);
        }
      }

      if (rank == size-1) fprintf(fd, "\n");

      fclose(fd);
    }
    PetscBarrier(NULL);
  }
} 


void VTKWriteFluidElementsParallel(char *VTKFile_Name, const PetscMPIInt size, const PetscMPIInt rank, const PetscInt Dim, 
  Nodes_Struct *nodes, Elements_Struct *elements, const PetscInt variable_number)
{  
  FILE *fd;
  PetscInt S_inter, S_bound, *GI;
  FLUID_ELEMENT_TYPE ET;

  if (variable_number == 1)      
    S_inter = elements->P_S_inter, S_bound = elements->P_S_bound, GI = nodes->GI_P, ET = elements->P_elem_type;
  else if (variable_number == 2) 
    S_inter = elements->U_S_inter, S_bound = elements->U_S_bound, GI = nodes->GI_U, ET = elements->U_elem_type;

  
  for (PetscInt r=0; r< size; ++r) 
  {
    if (rank == r) 
    {
      fd   = fopen(VTKFile_Name,"a");
      if(rank==0){ 
        if (ET != P1b)
          fprintf(fd, "CELLS %d %d\n", (elements->N_inter_glo+elements->N_bound_glo), 
            ((S_inter+1)*elements->N_inter_glo + (S_bound+1)*elements->N_bound_glo) );
        else  
          fprintf(fd, "CELLS %d %d\n", (elements->N_inter_glo*(Dim+1)+elements->N_bound_glo), 
            ((S_inter)*elements->N_inter_glo*(Dim+1) + (S_bound+1)*elements->N_bound_glo) );
      }
      if (variable_number == 2) 
        for(PetscInt i=0; i<elements->N_bound_loc; i++) {
          fprintf(fd, "%d", S_bound );
          for(PetscInt j=0; j<S_bound; j++)
            fprintf(fd, " %d", GI[elements->bound[i][j]] );
          fprintf(fd, "\n" );}
      else if (variable_number == 1) 
        for(PetscInt i=0; i<elements->N_bound_loc; i++) { 
          fprintf(fd, "%d", S_bound );
          for(PetscInt j=0; j<S_bound; j++)
            fprintf(fd, " %d", GI[nodes->LI_UtoP[elements->bound[i][j]]] );
          fprintf(fd, "\n" );}
      fclose(fd);
    }
    PetscBarrier(NULL);
  }

  
  for (PetscInt r=0; r< size; ++r) 
  {
    if (rank == r) 
     {
      fd   = fopen(VTKFile_Name,"a");

      if (ET == P1b)
      {
        if (Dim == 2)
          for(PetscInt i=0; i<elements->N_inter_loc; i++) 
          {
            fprintf(fd, "%d %d %d %d\n", S_inter-1, GI[elements->inter[i][0]], GI[elements->inter[i][1]], GI[elements->inter[i][3]]); 
            fprintf(fd, "%d %d %d %d\n", S_inter-1, GI[elements->inter[i][1]], GI[elements->inter[i][2]], GI[elements->inter[i][3]]); 
            fprintf(fd, "%d %d %d %d\n", S_inter-1, GI[elements->inter[i][2]], GI[elements->inter[i][0]], GI[elements->inter[i][3]]); 
          }
        else if (Dim == 3)
          for(PetscInt i=0; i<elements->N_inter_loc; i++) 
          {
            fprintf(fd, "%d %d %d %d %d\n", S_inter-1, GI[elements->inter[i][0]], GI[elements->inter[i][1]], GI[elements->inter[i][2]], GI[elements->inter[i][4]]); 
            fprintf(fd, "%d %d %d %d %d\n", S_inter-1, GI[elements->inter[i][1]], GI[elements->inter[i][2]], GI[elements->inter[i][3]], GI[elements->inter[i][4]]); 
            fprintf(fd, "%d %d %d %d %d\n", S_inter-1, GI[elements->inter[i][2]], GI[elements->inter[i][3]], GI[elements->inter[i][0]], GI[elements->inter[i][4]]); 
            fprintf(fd, "%d %d %d %d %d\n", S_inter-1, GI[elements->inter[i][3]], GI[elements->inter[i][0]], GI[elements->inter[i][1]], GI[elements->inter[i][4]]); 
          }
      }
      else if (variable_number == 2) 
        for(PetscInt i=0; i<elements->N_inter_loc; i++) 
        {
          fprintf(fd, "%d", S_inter);
          for(PetscInt j=0; j<S_inter; j++)
            fprintf(fd, " %d", GI[elements->inter[i][j]]);
          fprintf(fd, "\n" );
        }
      else if (variable_number == 1) 
        for(PetscInt i=0; i<elements->N_inter_loc; i++) 
        {
          fprintf(fd, "%d", S_inter);
          for(PetscInt j=0; j<S_inter; j++)
            fprintf(fd, " %d", GI[nodes->LI_UtoP[elements->inter[i][j]]]);
          fprintf(fd, "\n" );
        }

      if (rank == size-1) fprintf(fd, "\n" );

      fclose(fd);
    }
    PetscBarrier(NULL);
  }

  
  PetscInt VTK_inter, VTK_bound, Size_inter, Size_bound;
  if (ET!=P1b) VTKGetFluidElemInfo(&VTK_inter, &VTK_bound, &Size_inter, &Size_bound, Dim, ET);
  else         VTKGetFluidElemInfo(&VTK_inter, &VTK_bound, &Size_inter, &Size_bound, Dim, P1);

  if(rank==0)
  {
    PetscInt N_int, N_bou;
    if (ET!=P1b)  N_int = elements->N_inter_glo,          N_bou = elements->N_bound_glo;
    else          N_int = elements->N_inter_glo*(Dim+1),  N_bou = elements->N_bound_glo;

    fd   = fopen(VTKFile_Name,"a");
    fprintf(fd, "CELL_TYPES %d\n", N_int+N_bou); 
    for(PetscInt i=0; i<N_bou; i++)
      fprintf(fd, "%d\n", VTK_bound); 
    for(PetscInt i=0; i<N_int; i++)
      fprintf(fd, "%d\n", VTK_inter); 
    fprintf(fd, "\n" );
    fclose(fd);
  }
  PetscBarrier(NULL);
}

void VTKWritePointDataParallel_Ints(char *VTKFile_Name, const PetscMPIInt size, const PetscMPIInt rank, const char VarName[], 
  const PetscInt N_loc, const PetscInt N_glo, PetscInt *L, PetscInt *t,PetscBool AddHeader)
{
  FILE *fd;

  for (PetscInt r=0; r< size; ++r) 
  {
    if (rank == r)
    {
      fd = fopen(VTKFile_Name,"a");
      if (AddHeader && rank==0) fprintf(fd, "POINT_DATA %d\n", N_glo);
      if (rank ==0)             fprintf(fd, "SCALARS %s int 1\nLOOKUP_TABLE default\n", VarName);
      for(PetscInt i=0; i < N_loc; i++)
        if(t[i]<2) 
          fprintf(fd, "%d\n", L[i]);
      if (rank == size-1) fprintf(fd, "\n");
      fclose(fd);
    }
    PetscBarrier(NULL);
  }
}


PetscErrorCode VTKWriteSolidData(System_Struct *system, Solid_Struct *solid, const PetscInt time_step)
{
  PetscErrorCode ierr;
  char VTKFile_S_name[256] = "";

  const PetscMPIInt rank = system->rank;

  PetscMkdir(system->SYSOutputFolder);
  PetscSNPrintf(VTKFile_S_name,sizeof(VTKFile_S_name),"%s/S_%d.vtk",system->SYSOutputFolder, time_step);
  
  if (rank == 0)
  {
    VTKWriteIntro(VTKFile_S_name, "Solid");  
    VTKWriteNodes(VTKFile_S_name, solid->N_total_pts, solid->x, solid->y, solid->z);
    VTKWriteSolidElements(VTKFile_S_name, solid);
    VTKWritePointData_Vector(VTKFile_S_name, "F", solid->N_total_pts, solid->forces.Fx, solid->forces.Fy, solid->forces.Fz, PETSC_TRUE);
  }

  ierr = PetscBarrier(NULL); CHKERRQ(ierr);
  PetscPrintf(PETSC_COMM_WORLD,"  Solid file S_%d saved in VTK ASCII format\n", time_step);

  return 0;
}


void VTKWriteVectorSolution_U(char *VTKFile_Name, const PetscMPIInt size, const PetscMPIInt rank, const PetscInt Dim,
  Nodes_Struct *nodes, PetscInt **GI_V,const PetscScalar *V_values, const PetscInt i_low, const char *VectorName)
{
  FILE *fd;
  for (PetscInt r=0; r< size; ++r) 
  {
    if (rank == r)
    {
      fd = fopen(VTKFile_Name,"a");
      if (rank==0)
      {
        fprintf(fd, "POINT_DATA %d\n", nodes->N_U_total_glo);
        fprintf(fd, "VECTORS %s double\n", VectorName);
      }

      if (Dim == 2){
        for(PetscInt i=0; i<nodes->N_U_total_loc;++i) 
          if (nodes->t[i] < 2)      
            fprintf(fd, "%g %g 0.0\n", (double) V_values[GI_V[0][i]-i_low] , (double) V_values[GI_V[1][i]-i_low]); }
      else{
        for(PetscInt i=0; i<nodes->N_U_total_loc;++i) 
          if (nodes->t[i] < 2)      
            fprintf(fd, "%g %g %g\n", (double) V_values[GI_V[0][i]-i_low] , (double) V_values[GI_V[1][i]-i_low], (double) V_values[GI_V[2][i]-i_low]); }

      if (rank == size -1) fprintf(fd, "\n");
      fclose(fd);
    }
    PetscBarrier(NULL);
  }
}

void VTKWriteScalarSolution_P(char *VTKFile_Name, const PetscMPIInt size, const PetscMPIInt rank, const PetscInt Dim,
  Nodes_Struct *nodes, PetscInt *GI_V,const PetscScalar *V_values, const PetscInt i_low, const char *ScalarName)
{
  FILE *fd;
  for (PetscInt r=0; r< size; ++r) 
  {
    if (rank == r)
    {
      fd = fopen(VTKFile_Name,"a");
      if (rank==0)
      {
        fprintf(fd, "POINT_DATA %d\n", nodes->N_P_total_glo);
        fprintf(fd, "SCALARS %s double 1\n", ScalarName);
        fprintf(fd, "LOOKUP_TABLE default\n");
      }

     for(PetscInt i=0; i<nodes->N_P_total_loc;++i) 
        if (nodes->t[nodes->LI_PtoU[i]] < 2)  
          fprintf(fd, "%g\n", (double) V_values[GI_V[i]-i_low]);
         
      if (rank == size -1) fprintf(fd, "\n");
      fclose(fd);
    }
    PetscBarrier(NULL);
  }
}


PetscErrorCode VTKWriteFluidStokesVariables(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements,
  StokesVariables_Struct *stokesVariables, Vec *V, const PetscInt time_step)
{ 
  
  PetscErrorCode ierr;
  PetscInt V_size; 
  PetscScalar const *V_values; 
  PetscInt i_low; 
  char VTKFile_U_name[256] = "";
  char VTKFile_P_name[256] = "";

  const PetscMPIInt size = system->size;
  const PetscMPIInt rank = system->rank;
  const PetscInt    Dim  = system->Dim;

  VecGetArrayRead(*V, &V_values);   
  VecGetLocalSize(*V, &V_size);     
  VecGetOwnershipRange(*V, &i_low,NULL);
    
  PetscSNPrintf(VTKFile_U_name,sizeof(VTKFile_U_name),"%s/U_%d.vtk",system->SYSOutputFolder, time_step);
  VTKWriteIntroParallel(VTKFile_U_name, "Fluid-Velocity");  
  VTKWriteNodesParallel(VTKFile_U_name, size, rank, nodes, 2);
  VTKWriteFluidElementsParallel(VTKFile_U_name, size, rank, Dim, nodes, elements, 2);
  VTKWriteVectorSolution_U(VTKFile_U_name, size, rank, Dim, nodes, stokesVariables->GI_V, V_values, i_low, "U");
  PetscPrintf(PETSC_COMM_WORLD,"  Fluid vector file U_%d saved in VTK ASCII format\n", time_step);
  
  PetscSNPrintf(VTKFile_P_name,sizeof(VTKFile_P_name),"%s/P_%d.vtk",system->SYSOutputFolder, time_step);
  VTKWriteIntroParallel(VTKFile_P_name, "Fluid-Pressure");  
  VTKWriteNodesParallel(VTKFile_P_name, size, rank, nodes, 1);
  VTKWriteFluidElementsParallel(VTKFile_P_name, size, rank, Dim, nodes, elements, 1);
  VTKWriteScalarSolution_P(VTKFile_P_name, size, rank, Dim, nodes, stokesVariables->GI_V[stokesVariables->Nv_U], V_values, i_low, "P");
  PetscPrintf(PETSC_COMM_WORLD,"  Fluid scalar file P_%d saved in VTK ASCII format\n", time_step);

  VecRestoreArrayRead(*V,&V_values);
  
  ierr = PetscBarrier(NULL); CHKERRQ(ierr);
  return 0;
}


PetscErrorCode VTKWriteFluidMesh(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements)
{
  const PetscMPIInt size = system->size;
  const PetscMPIInt rank = system->rank;
  const PetscInt    Dim  = system->Dim;
  char VTKFile_Name[256] = "";

  if (size > 1)
  {
    
    PetscSNPrintf(VTKFile_Name,sizeof(VTKFile_Name),"%s/FluidMesh_U_%d.vtk",system->SYSOutputFolder, rank);
    VTKWriteIntro(VTKFile_Name, "FluidMesh_U");  
    VTKWriteNodes(VTKFile_Name, nodes->N_U_total_loc, nodes->x, nodes->y, nodes->z);
    VTKWriteFluidElements(VTKFile_Name, Dim, nodes, elements, 2);
    VTKWritePointData_Ints(VTKFile_Name, "t",  nodes->N_U_total_loc, nodes->t,    PETSC_TRUE);
    VTKWritePointData_Ints(VTKFile_Name, "GI", nodes->N_U_total_loc, nodes->GI_U, PETSC_FALSE);
    
    PetscSNPrintf(VTKFile_Name,sizeof(VTKFile_Name),"%s/FluidMesh_P_%d.vtk",system->SYSOutputFolder, rank);
    VTKWriteIntro(VTKFile_Name, "FluidMesh_P");  
    VTKWriteNodesP(VTKFile_Name, nodes->N_P_total_loc, nodes->x, nodes->y, nodes->z, nodes->LI_PtoU);
    VTKWriteFluidElements(VTKFile_Name, Dim, nodes, elements, 1);
    VTKWritePointDataP_Ints(VTKFile_Name, "t",  nodes->N_P_total_loc, nodes->t, nodes->LI_PtoU, PETSC_TRUE);
    VTKWritePointData_Ints(VTKFile_Name, "GI", nodes->N_P_total_loc, nodes->GI_P, PETSC_FALSE);
  }

  PetscSNPrintf(VTKFile_Name,sizeof(VTKFile_Name),"%s/FluidMesh_U.vtk",system->SYSOutputFolder);
  VTKWriteIntroParallel(VTKFile_Name, "Fluid-Velocity");  
  VTKWriteNodesParallel(VTKFile_Name, size, rank, nodes, 2);
  VTKWriteFluidElementsParallel(VTKFile_Name, size, rank, Dim, nodes, elements, 2);
  VTKWritePointDataParallel_Ints(VTKFile_Name, size, rank,  "t", nodes->N_U_total_loc, nodes->N_U_total_glo, nodes->t, nodes->t, PETSC_TRUE);
  VTKWritePointDataParallel_Ints(VTKFile_Name, size, rank, "GI", nodes->N_U_total_loc, nodes->N_U_total_glo, nodes->GI_U, nodes->t, PETSC_FALSE);

  PetscSNPrintf(VTKFile_Name,sizeof(VTKFile_Name),"%s/FluidMesh_P.vtk",system->SYSOutputFolder);
  VTKWriteIntroParallel(VTKFile_Name, "Fluid-Pressure");  
  VTKWriteNodesParallel(VTKFile_Name, size, rank, nodes, 1);
  VTKWriteFluidElementsParallel(VTKFile_Name, size, rank, Dim, nodes, elements, 1);
  return 0;
}


PetscErrorCode VTKCreateSolidMesh(System_Struct *system, Solid_Struct *solid)
{
  PetscErrorCode ierr;
  PetscLogDouble Clock1, Clock2;
  char SDVTK_File[256] = "";
  char line[256]="";
  FILE *fd;

  PetscInt VTKPOINTS, VTKLINES, VTKTRIANGLES;
  VTKGetSolidElemInfo(&VTKPOINTS, POINT);
  VTKGetSolidElemInfo(&VTKLINES, LINE);
  VTKGetSolidElemInfo(&VTKTRIANGLES, TRIANGLE);

  
  PetscBool flg_POINTS      = PETSC_FALSE;
  PetscBool flg_CELLS       = PETSC_FALSE;
  PetscBool flg_CELL_TYPES  = PETSC_FALSE;
  PetscBool flg_CELL_DATA   = PETSC_FALSE;

  PetscInt *cell_groups;
  PetscInt *bound_cells, *bound_cell_type, *bound_groups;
  PetscInt cells_total, cells_types_total, cell_data_total=-1;
  PetscInt k;
  char **data;
  char data2[256];

  PetscInt *ln0, *ln1, *tr0, *tr1, *tr2;

  PetscTime(&Clock1);

  PetscPrintf(PETSC_COMM_WORLD,"Reading VTK '.vtk' files on all ranks\n");
  PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------\n"); 

  PetscSNPrintf(SDVTK_File,sizeof(SDVTK_File),"%s.vtk",system->SDMesh_File);

  PetscPrintf(PETSC_COMM_WORLD,"  VTK File: %s\n", SDVTK_File);
  
  if (system->rank == 0)
  {
    fd   = fopen(SDVTK_File,"r");

    if(fd == NULL) 
    { PetscErrorPrintf("VTK '.vtk' file %s is not found.\n", SDVTK_File);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }
    
    while (PETSC_TRUE)
    {
      ReadFileLine(fd, 256,line); 
      PetscStrncmp(line, "POINTS",           6, &flg_POINTS);

      if(flg_POINTS)
      { 
        PetscStrToArray(line,' ',&k, &data);
        ierr = PetscOptionsStringToInt(data[1], &solid->N_total_pts);  CHKERRQ(ierr); 
        PetscPrintf(PETSC_COMM_WORLD,"  Total number of points: %d\n",solid->N_total_pts);
        PetscStrToArrayDestroy(k, data);

        SolidMallocPointArrays(solid);
        
        for(PetscInt i=0; i < solid->N_total_pts; i++)
        {
          ReadFileLine(fd, 256,line);                                    
          PetscStrToArray(line,' ',&k, &data);
          ierr = PetscOptionsStringToScalar(data[0], &solid->x[i]);    CHKERRQ(ierr); 
          ierr = PetscOptionsStringToScalar(data[1], &solid->y[i]);    CHKERRQ(ierr); 
          RemoveNewLine(data[2], data2); 
          ierr = PetscOptionsStringToScalar(data2, &solid->z[i]);      CHKERRQ(ierr);
          PetscStrToArrayDestroy(k, data);
        }
        break;
      }

      if (line[0] == 0)  break; 
    }
    
    ierr = fclose(fd);
    if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

    if (system->size > 1)
    {
      ierr = MPI_Bcast(&solid->N_total_pts, 1, MPIU_INT, 0, PETSC_COMM_WORLD );       CHKERRMPI(ierr); 
      ierr = MPI_Bcast(solid->x, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
      ierr = MPI_Bcast(solid->y, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
      ierr = MPI_Bcast(solid->z, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    }
  }
  else
  {
    ierr = MPI_Bcast(&solid->N_total_pts, 1, MPIU_INT, 0, PETSC_COMM_WORLD );       CHKERRMPI(ierr); 
    SolidMallocPointArrays(solid);
    ierr = MPI_Bcast(solid->x, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
    ierr = MPI_Bcast(solid->y, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
    ierr = MPI_Bcast(solid->z, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
  }

  PetscBarrier(NULL);
  PetscPrintf(PETSC_COMM_WORLD,"  Done reading solid POINTS section\n");
  
  if (system->rank == 0)
  {
    fd = fopen(SDVTK_File,"r");

    if(fd == NULL) 
    { PetscErrorPrintf("VTK '.vtk' file %s is not found.\n", SDVTK_File);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

    while (PETSC_TRUE)
    {
      ReadFileLine(fd, 256,line);  
      PetscStrncmp(line, "CELL_DATA",  9, &flg_CELL_DATA);

      if (flg_CELL_DATA)
      {
        PetscStrToArray(line,' ',&k, &data);
        RemoveNewLine(data[1], data2); 
        ierr = PetscOptionsStringToInt(data2, &cell_data_total);      CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD,"  Total number of cell: %d\n",cell_data_total);
        PetscStrToArrayDestroy(k, data);
        ReadFileLine(fd, 256,line);     
        ReadFileLine(fd, 256,line);     

        ierr = PetscMalloc1(cell_data_total, &cell_groups); CHKERRQ(ierr);

        
        for(PetscInt i=0; i < cell_data_total; i++)
        {
          ReadFileLine(fd, 256,line);                                    
          PetscStrToArray(line,' ',&k, &data);
          RemoveNewLine(data[0], data2); 
          ierr = PetscOptionsStringToInt(data2, &cell_groups[i]); CHKERRQ(ierr);
          PetscStrToArrayDestroy(k, data);
        }
        break;
      }

      if (line[0] == 0)  break; 
    }
    ierr = fclose(fd);
    if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

    if (cell_data_total==-1) 
    {
      solid->N_bound = 0;
      solid->N_Gr = 0;
    }
    else
    {
      solid->N_bound = 0;
      for(PetscInt i=0; i < cell_data_total; i++)
        if (cell_groups[i] != -1) solid->N_bound++;

      solid->N_Gr = 0;
      if (solid->N_bound != 0)
      {
        PetscInt MIN_GR = 100000, MAX_GR = -10000;
        for (PetscInt i=0; i< cell_data_total;i++){
          const PetscInt bg = cell_groups[i];
          if (bg > MAX_GR && bg != -1) MAX_GR = bg;
          if (bg < MIN_GR && bg != -1) MIN_GR = bg;
        }
        
        solid->N_Gr = MAX_GR - MIN_GR + 1;
      }
      PetscPrintf(PETSC_COMM_WORLD,"  The number of boundary groups is %d\n",solid->N_Gr);
    }

    if (system->size>1){
      ierr = MPI_Bcast(&solid->N_Gr, 1, MPIU_INT, 0, PETSC_COMM_WORLD );       CHKERRMPI(ierr);}

    PetscMalloc3(solid->N_bound, &bound_cells, solid->N_bound, &bound_cell_type, solid->N_bound, &bound_groups);
  }
  else
  {
    ierr = MPI_Bcast(&solid->N_Gr, 1, MPIU_INT, 0, PETSC_COMM_WORLD );       CHKERRMPI(ierr);
  }

  SolidMallocBoundaryGroupArrays(solid);
  PetscBarrier(NULL);
  PetscPrintf(PETSC_COMM_WORLD,"  Done reading solid CELL_DATA section. Number of boundary elements = %d\n", solid->N_bound);

  
  if (system->rank == 0)
  {
    PetscInt cell_t, b=0;

    solid->N_total_lns=0; 
    solid->N_total_trs=0; 

    fd   = fopen(SDVTK_File,"r");

    if(fd == NULL) 
    { PetscErrorPrintf("VTK '.vtk' file %s is not found.\n", SDVTK_File);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

    
    while (PETSC_TRUE)
    {
      
      ReadFileLine(fd, 256,line); 
      PetscStrncmp(line, "CELL_TYPES",          10, &flg_CELL_TYPES);

      if(flg_CELL_TYPES)
      {
        PetscStrToArray(line,' ',&k, &data);
        RemoveNewLine(data[1], data2); 
        ierr = PetscOptionsStringToInt(data2, &cells_types_total);      CHKERRQ(ierr);
        PetscStrToArrayDestroy(k, data);

        if(cell_data_total!=-1)
          if (cells_types_total != cell_data_total)    
            { PetscErrorPrintf("Error cells_types_total:%d != cell_data_total:%d\n", cells_types_total, cell_data_total);
              SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in reading CELL_TYPES section"); }

        for(PetscInt i=0; i < cells_types_total; i++)
        {
          ReadFileLine(fd, 256,line);                                    
          PetscStrToArray(line,' ',&k, &data);
          RemoveNewLine(data[0], data2); 
          ierr = PetscOptionsStringToInt(data2, &cell_t);      CHKERRQ(ierr);
          if (cell_t == VTKLINES)           solid->N_total_lns++;
          else if (cell_t == VTKTRIANGLES)  solid->N_total_trs++;
          PetscStrToArrayDestroy(k, data);

          
          if (cell_data_total!=-1 && cell_groups[i] != -1)
          {
            bound_cell_type[b] = cell_t;
            bound_cells[b] = i;
            bound_groups[b] = cell_groups[i] - 1;
            b++;
          }

        }
        break;

      }
      if (line[0] == 0)  break; 
    }
    
    ierr = fclose(fd);
    if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

    
    solid->N_bound_pts=0, solid->N_bound_lns=0, solid->N_bound_trs=0;
    for (PetscInt b=0; b<solid->N_bound; b++)
    {
      const PetscInt g = bound_groups[b]; 

      if (bound_cell_type[b] == VTKPOINTS)    {
        solid->N_bound_pts++;
        if (solid->Gr_elem_type[g] == NO_SET) solid->Gr_elem_type[g] = POINT;
        else if (solid->Gr_elem_type[g] != POINT){
          PetscErrorPrintf("Error! Each solid boundary can only have one element type. %s != %s\n", 
            SOLID_ELEMENT_TYPE_CHAR[solid->Gr_elem_type[g]], SOLID_ELEMENT_TYPE_CHAR[POINT]);
          SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Error in reading CELL_TYPES section for boundary");}
      }
      else if (bound_cell_type[b] == VTKLINES)     {
        solid->N_bound_lns++;
        if (solid->Gr_elem_type[g] == NO_SET) solid->Gr_elem_type[g] = LINE;
        else if (solid->Gr_elem_type[g] != LINE){
          PetscErrorPrintf("Error! Each solid boundary can only have one element type. %s != %s\n", 
            SOLID_ELEMENT_TYPE_CHAR[solid->Gr_elem_type[g]], SOLID_ELEMENT_TYPE_CHAR[LINE]);
          SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Error in reading CELL_TYPES section for boundary");}
      }
      else if (bound_cell_type[b] == VTKTRIANGLES) {
        solid->N_bound_trs++; 
        if (solid->Gr_elem_type[g] == NO_SET) solid->Gr_elem_type[g] = TRIANGLE;
        else if (solid->Gr_elem_type[g] != TRIANGLE){
          PetscErrorPrintf("Error! Each solid boundary can only have one element type. %s != %s\n", 
            SOLID_ELEMENT_TYPE_CHAR[solid->Gr_elem_type[g]], SOLID_ELEMENT_TYPE_CHAR[TRIANGLE]);
          SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Error in reading CELL_TYPES section for boundary");}
      }
    }
      
    if (system->size > 1)
    {
      ierr = MPI_Bcast(&solid->N_total_lns, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(&solid->N_total_trs, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(&solid->N_bound,     1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(&solid->N_bound_pts, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(&solid->N_bound_lns, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(&solid->N_bound_trs, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    }
  }
  else
  {
    ierr = MPI_Bcast(&solid->N_total_lns, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    ierr = MPI_Bcast(&solid->N_total_trs, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    ierr = MPI_Bcast(&solid->N_bound,     1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    ierr = MPI_Bcast(&solid->N_bound_pts, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    ierr = MPI_Bcast(&solid->N_bound_lns, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    ierr = MPI_Bcast(&solid->N_bound_trs, 1, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
  }

  SolidMallocCellArrays(solid);
  SolidMallocBoundaryCellArrays(solid);
  PetscBarrier(NULL);
  PetscPrintf(PETSC_COMM_WORLD,"  Done reading solid CELL_TYPES section\n");
  PetscPrintf(PETSC_COMM_WORLD,"    Number of lines = %d\n    Number of triangles = %d\n    Number of boundary cells = %d\n", 
    solid->N_total_lns, solid->N_total_trs, solid->N_bound);
  PetscPrintf(PETSC_COMM_WORLD,"    Number of boundary points = %d\n    Number of boundary lines = %d\n    Number of boundary triangles = %d\n", 
    solid->N_bound_pts, solid->N_bound_lns, solid->N_bound_trs);

  if (system->rank == 0)
  {
    PetscInt cell_size, l=0, t=0, bp=0, bl=0, bt=0;

    fd   = fopen(SDVTK_File,"r");

    if(fd == NULL) 
    { PetscErrorPrintf("VTK '.vtk' file %s is not found.\n", SDVTK_File);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

    
    while (PETSC_TRUE)
    {
      ReadFileLine(fd, 256,line);  
      PetscStrncmp(line, "CELLS",  5, &flg_CELLS);

      if(flg_CELLS)
      {
        PetscStrToArray(line,' ',&k, &data);
        ierr = PetscOptionsStringToInt(data[1], &cells_total);      CHKERRQ(ierr);
        PetscStrToArrayDestroy(k, data);
        
        if (cells_total != cells_types_total)    
          { PetscErrorPrintf("Error cells_total:%d != cells_types_total:%d\n", cells_total, cells_types_total);
            SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error in reading CELLS section."); }

        
        for(PetscInt i=0; i < cells_total; i++)
        {
          ReadFileLine(fd, 256,line);                                    
          PetscStrToArray(line,' ',&k, &data);
          ierr = PetscOptionsStringToInt(data[0], &cell_size); CHKERRQ(ierr);
          if (cell_size == 2) { 
            ierr = PetscOptionsStringToInt(data[1], &solid->lns[l][0]); CHKERRQ(ierr);
            RemoveNewLine(data[2], data2); 
            ierr = PetscOptionsStringToInt(data2, &solid->lns[l][1]);   CHKERRQ(ierr);
            l++; }
          else if (cell_size == 3) { 
            ierr = PetscOptionsStringToInt(data[1], &solid->trs[t][0]); CHKERRQ(ierr);
            ierr = PetscOptionsStringToInt(data[2], &solid->trs[t][1]); CHKERRQ(ierr);
            RemoveNewLine(data[3], data2); 
            ierr = PetscOptionsStringToInt(data2, &solid->trs[t][2]);   CHKERRQ(ierr); 
            t++;}

          
          if (cell_data_total!=-1 && cell_groups[i] != -1){
            if (cell_size == 1){
              RemoveNewLine(data[1], data2); 
              ierr = PetscOptionsStringToInt(data2, &solid->bound_pts[bp]); CHKERRQ(ierr); 
              solid->bound_group_pts[bp] = cell_groups[i]-1;   bp++; }
            else if (cell_size == 2){
              solid->bound_lns[bl] = l-1; solid->bound_group_lns[bl] = cell_groups[i]-1; bl++; }
            else if (cell_size == 3){
              solid->bound_trs[bt] = t-1; solid->bound_group_trs[bt] = cell_groups[i]-1; bt++; }
          }

          PetscStrToArrayDestroy(k, data);
        }
        break;
      }
      if (line[0] == 0)  break; 
    }
    
    ierr = fclose(fd);
    if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

    
    if (system->size > 1)
    { 
      if (solid->N_total_lns>0)
      {
        PetscMalloc2(solid->N_total_lns, &ln0, solid->N_total_lns, &ln1);
        for (PetscInt i=0; i<solid->N_total_lns; i++)
          ln0[i] = solid->lns[i][0], ln1[i] = solid->lns[i][1];
        ierr = MPI_Bcast(ln0, solid->N_total_lns, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
        ierr = MPI_Bcast(ln1, solid->N_total_lns, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
        PetscFree2(ln0, ln1);
      }

      if (solid->N_total_trs>0)
      {
        PetscMalloc3(solid->N_total_trs, &tr0, solid->N_total_trs, &tr1, solid->N_total_trs, &tr2);
        for (PetscInt i=0; i<solid->N_total_trs; i++)
          tr0[i] = solid->trs[i][0], tr1[i] = solid->trs[i][1], tr2[i] = solid->trs[i][2];
        ierr = MPI_Bcast(tr0, solid->N_total_trs, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
        ierr = MPI_Bcast(tr1, solid->N_total_trs, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
        ierr = MPI_Bcast(tr2, solid->N_total_trs, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
        PetscFree3(tr0, tr1, tr2);
      }

      if (solid->N_bound > 0)
      {
        ierr = MPI_Bcast(solid->bound_pts,       solid->N_bound_pts, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
        ierr = MPI_Bcast(solid->bound_group_pts, solid->N_bound_pts, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
        ierr = MPI_Bcast(solid->bound_lns,       solid->N_bound_lns, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
        ierr = MPI_Bcast(solid->bound_group_lns, solid->N_bound_lns, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
        ierr = MPI_Bcast(solid->bound_trs,       solid->N_bound_trs, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
        ierr = MPI_Bcast(solid->bound_group_trs, solid->N_bound_trs, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      }
    }
  }
  else
  {
    
    if (solid->N_total_lns>0)
    {
      PetscMalloc2(solid->N_total_lns, &ln0, solid->N_total_lns, &ln1);
      ierr = MPI_Bcast(ln0, solid->N_total_lns, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(ln1, solid->N_total_lns, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      for (PetscInt i=0; i<solid->N_total_lns; i++)
        solid->lns[i][0] = ln0[i], solid->lns[i][1] = ln1[i];
      ierr = PetscFree2(ln0, ln1); CHKERRQ(ierr);
    }

    if (solid->N_total_trs>0)
    {
      PetscMalloc3(solid->N_total_trs, &tr0, solid->N_total_trs, &tr1, solid->N_total_trs, &tr2);
      ierr = MPI_Bcast(tr0, solid->N_total_trs, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(tr1, solid->N_total_trs, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(tr2, solid->N_total_trs, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      for (PetscInt i=0; i<solid->N_total_trs; i++)
        solid->trs[i][0] = tr0[i], solid->trs[i][1] = tr1[i], solid->trs[i][2] = tr2[i];
      PetscFree3(tr0, tr1, tr2); 
    }

    if (solid->N_bound > 0)
    {
      ierr = MPI_Bcast(solid->bound_pts,       solid->N_bound_pts, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(solid->bound_group_pts, solid->N_bound_pts, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(solid->bound_lns,       solid->N_bound_lns, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(solid->bound_group_lns, solid->N_bound_lns, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(solid->bound_trs,       solid->N_bound_trs, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
      ierr = MPI_Bcast(solid->bound_group_trs, solid->N_bound_trs, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    }
  }

  if(system->rank==0){
    if(solid->N_bound >0 ) {ierr = PetscFree(cell_groups); CHKERRQ(ierr);}
    ierr = PetscFree3(bound_cells, bound_cell_type, bound_groups);CHKERRQ(ierr);
  }

  PetscPrintf(PETSC_COMM_WORLD,"  Done reading the solid points VTK file\n"); 
  ierr = PetscBarrier(NULL); CHKERRQ(ierr);
  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"  Computational time: %lf seconds\n\n", (Clock2 - Clock1) ); 

  return ierr;
}



PetscErrorCode ReadVTKSolidFileAndUpdatePoints(System_Struct *system, Solid_Struct *solid)
{
  
  PetscErrorCode ierr;
  PetscLogDouble Clock1, Clock2;
  char SDVTK_File[256] = "";
  char line[256]="";
  FILE *fd;
  PetscInt N_total_pts_Reload = 0;

  PetscBool flg_POINTS      = PETSC_FALSE;

  PetscInt k;
  char **data;
  char data2[256];

  PetscTime(&Clock1);

  PetscPrintf(PETSC_COMM_WORLD,"Reading VTK '.vtk' files on all ranks\n");

  PetscSNPrintf(SDVTK_File,sizeof(SDVTK_File),"%s.vtk",system->SDReload_File);

  PetscPrintf(PETSC_COMM_WORLD,"  VTK File: %s\n", SDVTK_File);

  if (system->rank == 0)
  {
    fd   = fopen(SDVTK_File,"r");

    if(fd == NULL) 
    { PetscErrorPrintf("VTK '.vtk' file %s is not found.\n", SDVTK_File);
      SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

    while (PETSC_TRUE)
    {
      ReadFileLine(fd, 256,line); 
      PetscStrncmp(line, "POINTS",           6, &flg_POINTS);

      if(flg_POINTS)
      { 
        PetscStrToArray(line,' ',&k, &data);
        ierr = PetscOptionsStringToInt(data[1], &N_total_pts_Reload);  CHKERRQ(ierr); 
        PetscStrToArrayDestroy(k, data);

        if(N_total_pts_Reload != solid->N_total_pts) { 
          PetscErrorPrintf("The number of points in the reload file do not match the original number of points. i.e. %d != %d .\n", N_total_pts_Reload, solid->N_total_pts);
          SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! Wrong reload file"); }

        for(PetscInt i=0; i < N_total_pts_Reload; i++)
        {
          ReadFileLine(fd, 256,line);                                    
          PetscStrToArray(line,' ',&k, &data);
          ierr = PetscOptionsStringToScalar(data[0], &solid->x[i]);    CHKERRQ(ierr); 
          ierr = PetscOptionsStringToScalar(data[1], &solid->y[i]);    CHKERRQ(ierr); 
          RemoveNewLine(data[2], data2); 
          ierr = PetscOptionsStringToScalar(data2, &solid->z[i]);      CHKERRQ(ierr);
          PetscStrToArrayDestroy(k, data);
        }
        break;
      }

      if (line[0] == 0)  break; 
    }
    
    ierr = fclose(fd);
    if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

    if (system->size > 1)
    {
      ierr = MPI_Bcast(solid->x, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); PetscBarrier(NULL);
      ierr = MPI_Bcast(solid->y, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); PetscBarrier(NULL);
      ierr = MPI_Bcast(solid->z, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    }
  }
  else
  {
    ierr = MPI_Bcast(solid->x, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); PetscBarrier(NULL);
    ierr = MPI_Bcast(solid->y, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); PetscBarrier(NULL);
    ierr = MPI_Bcast(solid->z, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
  }

  PetscBarrier(NULL);
  PetscPrintf(PETSC_COMM_WORLD,"  Done reloading solid POINTS\n");
  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"  Computational time: %lf seconds\n\n", (Clock2 - Clock1) ); 

  return 0;
}

#endif