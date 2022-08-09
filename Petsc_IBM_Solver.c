static char help[] = "Petsc Stokes equation solver.\n\n";

// Petsc Headers:
//====================================================
#include <petscksp.h>
#include <petsclog.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscdmplex.h>

#include "Functions/Misc/Misc.c"
#include "Functions/Math/MathFunctions.c"

// Global Variables:
//====================================================
#include "Global_Files/Global_Enum.h"
#include "Global_Files/Global_Variables.h"
#include "Objects/System/SystemStruct.h"
#include "Objects/FluidMesh/NodesStruct.h"
#include "Objects/FluidMesh/ElementsStruct.h"
#include "Objects/StokesVariables/StokesVariablesStruct.h"
#include "Objects/Boxes/BoxesStruct.h"
#include "Objects/Solid/SolidStruct.h"

// Functions:
//====================================================
#include "Functions/FEM/FEMFunctions.h"
#include "Objects/Solid/SolidFunctions.c"
#include "Functions/ReadWrite/ReadWriteFunctions.h"
#include "Functions/Mesh/MeshFunctions.c"
#include "Functions/LinearSystem/LinearSystemFunctions.h"
#include "Functions/KSPSystem/KSPSystemFunctions.c"
#include "Objects/StokesVariables/StokesVariablesFunctions.c"
#include "Objects/Boxes/BoxesFunctions.c"
#include "Objects/Solid/SolidReadBoundary.c"
#include "Objects/Solid/Forces/ForcesFunctions.c"
#include "Functions/IB/IBFunctions.c"

// Main Function:
//====================================================
int main(int argc, char **args) 
{
  // Global Clock
  PetscLogDouble ClockStart, ClockEnd; 

  // Petsc Variables
  PetscTime(&ClockStart);
  #include "Header_Files/Initial_Variables.h"

  //====================================================
  // Program Starting:
  //====================================================

  // Initialize
  ierr = PetscInitialize(&argc,&args,(char*)0,help); if (ierr) return ierr;
  ierr = MPI_Comm_size(PETSC_COMM_WORLD, &system.size);CHKERRMPI(ierr);
  ierr = MPI_Comm_rank(PETSC_COMM_WORLD, &system.rank);CHKERRMPI(ierr);

  stdout_viewer = PETSC_VIEWER_STDOUT_WORLD;

  PetscPrintf(PETSC_COMM_WORLD,"\n\nPETSc IBM Solver\n"); 
  PetscPrintf(PETSC_COMM_WORLD,"==============================\n\n"); 

  #include "Header_Files/Inline_Input.h"
  #include "Header_Files/InputFile_Parameters.h"

  CreateFluidDomain(&system, &nodes, &elements);

  CreateSolidDomain(&system, &solid);
  SolidReadBoundary(&system, &solid);
  SolidInitializeProperties(&system, &solid);

  BoxesCreateDomain(&system, &nodes, &elements, &boxes);
  SearchForBoxWithPoints(&system, &nodes, &elements, &boxes, &solid);

  StokesVariablesCreate(&system, &nodes, &elements, &stokesVariables);

  AssembleStokesFullLS(&system, &nodes, &elements, &stokesVariables, &stokesFullLS);
  if(system.FLSolver_Method == MONO || system.FLSolver_Method == MONO_PRE)
    ConvertStokesFullLSToReducedMonoLS(&system, &nodes, &elements, &stokesVariables, &stokesFullLS, &stokesMonoLS);
  if(system.FLSolver_Method == MONO_PRE)
    CreateMonolithicPreconditionarMatrix(&system, &nodes, &elements, &stokesVariables, &stokesFullLS, &stokesMonoLS);

  CreateKSPForStokesLS(&system, &stokesFullLS, &stokesMonoLS, &ksp, &pc);

  SolveKSPForStokesLS(&system, &stokesFullLS, &stokesMonoLS, &ksp, &pc);

  IBTimeStepping(&system, &nodes, &elements, &stokesVariables, 
    &stokesFullLS, &stokesMonoLS, &boxes, &solid, &active_forces, &ksp, &pc);

  PetscPrintf(PETSC_COMM_WORLD,"\nCleaning memory\n"); 
  PetscPrintf(PETSC_COMM_WORLD,"---------------------------------------------\n"); 

  ReducedMonoLinearSystem_Stokes_Free(&system, &stokesMonoLS);
  FullLinearSystem_Stokes_Free(&system, &stokesFullLS);
  KSPDestroy(&ksp);
  StokesVariablesFree(&stokesVariables);
  BoxesFree(&boxes);
  ElementsFree(&elements); 
  NodesFree(&nodes, system.size);
  ForcesFree(&solid);
  SolidFree(&solid);

  PetscTime(&ClockEnd);
  PetscPrintf(PETSC_COMM_WORLD,"\nSimulation Finished\nSimulation time: %lf seconds\n", (ClockEnd - ClockStart));
  ierr = PetscFinalize();
  return ierr;
}
