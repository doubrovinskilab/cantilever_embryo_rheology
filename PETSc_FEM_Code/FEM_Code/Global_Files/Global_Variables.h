
PetscViewer  stdout_viewer; 

PetscBool    print_mat_A       = PETSC_FALSE; // Global A (Full LS)
PetscBool    print_vec_b       = PETSC_FALSE; // Global b (Full LS)

PetscBool    print_mat_Ar      = PETSC_FALSE; // Global Ar (Mono Reduced LS)
PetscBool    print_mat_Pr      = PETSC_FALSE; // Global Pr (Mono Reduced Preconditioner LS)
PetscBool    print_vec_br      = PETSC_FALSE; // Global br (Mono Reduced LS)
PetscBool    print_vec_b_dir   = PETSC_FALSE; // Global b_dir 
PetscBool    check_symmetry_Ar = PETSC_FALSE; // Check symmetry of the Mono reduced matrix

PetscBool    print_dir_V_U     = PETSC_FALSE;
PetscBool    print_dir_V_P     = PETSC_FALSE;
PetscBool    print_non_dir_V_U = PETSC_FALSE;
PetscBool    print_non_dir_V_P = PETSC_FALSE;

PetscBool    print_nodes            = PETSC_FALSE;
PetscBool    print_elements         = PETSC_FALSE;
PetscBool    print_stokesVariables  = PETSC_FALSE;
PetscBool    print_solid            = PETSC_FALSE;

PetscBool    write_solid_fe         = PETSC_FALSE; // Add elastic force to the vtk of solid

PetscBool    write_fluid_mesh       = PETSC_FALSE;
PetscBool    write_fluid_info       = PETSC_FALSE;
PetscInt     write_fluid_N_dtime    = 0; // Based on dtime in system 

const PetscScalar  EPS = 1e-8; // Epsilon
const PetscScalar   PI = 3.14159265358979323846; // PI
