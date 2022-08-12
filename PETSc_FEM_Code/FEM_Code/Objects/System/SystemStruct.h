
typedef struct {

  PetscMPIInt size, rank;                       // MPI parameters
  PetscInt    Dim;                              // Domain Dimension

  // Time parameters:
  //===============================
  PetscScalar         time;                     // Real simulation time
  PetscScalar         dtime_save;
  PetscScalar         end_time;
  PetscScalar         dt;                       // Either FLImp_dt or SD_dt

  // System Options:
  //===============================
  SYSTEM_TYPE         SYSType;                  // STOKES or IBM
  char                SYSOutputFolder[256];     // Folder to save
  char                SYSInputFile[256];        // File to read input from

  // Fluid Options:
  //===============================

  // Solver Options:
  FLUID_SOLVER_TYPE   FLSolver_Type;             // { Implicit, Explicit, Explicit_Lumped, Implicit_Explicit_Lumped }
  FLUID_SOLVER_METHOD FLSolver_Method;           // { Direct, Nest, Mono, Mono_Pre}
  PetscBool           FLSolver_AddHHCorr;        // Solve Helmholtz-Hodge Correction
  PetscBool           FLSolver_DivVCorr;         // Solve the Div.U=0 at the solid surface
  PetscBool           FLSolver_ForceAddVolCons;  // Adds a volume conservation force

  // Implicit Parameters
  PetscScalar         FLImp_dt;                  // dt to be used when solving the implict equation
  PetscScalar         FLImp_eps;                 // eps to be used when solving the implict penalty equation (if eps=0 ==> solve regular equation)

  // Explicit Parameters
  PetscScalar         FLExp_dt;                  // dt to be used when solving the explicit fluid solution      
  PetscInt            FLExp_N;                   // N number of iterations when solving the explicit fluid solution
  PetscScalar         FLExp_eps;                 // eps to be used when solving the explicit penalty equation
  PetscScalar         FLExp_rho;                // rho (Default = 1)

  // Parameters for Solving Fluid implicitly/explicitly
  PetscScalar         FLImpExp_rTol;             // Tolarence to be used when switching between Imp or Explicit Schemes
  PetscBool           FLImpExp_Bool;             // A boolean to tell the code if implicit method was utilized in the current time step. True=> Imp, False => Exp

  // Mesh Options:
  MESH_TYPE           FLMesh_Type;               // { vtk, gmsh}
  char                FLMesh_File[256];          // Mesh file name 

  // Solid Options:
  //===============================
  // Mesh Options:
  PetscScalar               SD_dt;                // dt to be used when solving the linear ealsicity only
  MESH_TYPE           SDMesh_Type;                // { vtk, gmsh}
  char                SDMesh_File[256];           // Mesh file name 
  PetscBool           SDMesh_ConvTrs2Lns;        // A boolean to convert triangle cells to line cells

  PetscBool           SDDivFree_Bool;             // A boolean to tell the code to convert the solid velocity field to a divergent field

  PetscBool           SDReload;                   // An option to reload a data file
  char                SDReload_File[256];         // A solid file that is used to reload a configuration and work from it
  PetscScalar         SDReload_time ;             // Reloaded time
  PetscInt            SDReload_saved_step;        // Reloaded saved step


} System_Struct; // Linear System

#include "SystemStructFunctions.c"
