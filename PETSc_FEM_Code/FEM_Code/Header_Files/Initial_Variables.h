
PetscErrorCode  ierr;
PetscBool       flg;

System_Struct               system;
Nodes_Struct                nodes;
Elements_Struct             elements;
StokesVariables_Struct      stokesVariables;
Solid_Struct                solid;
Boxes_Struct                boxes;
ForcesActive_Struct         active_forces;

FullLinearSystem_Stokes_Struct stokesFullLS;
ReducedMonoLinearSystem_Stokes_Struct stokesMonoLS;

KSP ksp;
PC  pc;

// Initialize Structs
SystemInitialize(&system);
NodesInitialize(&nodes);
ElementsInitialize(&elements);
StokesVariablesInitialize(&stokesVariables);
SolidInitialize(&solid);
BoxesInitialize(&boxes);