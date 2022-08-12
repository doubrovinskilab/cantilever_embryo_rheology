
  typedef struct {
    PetscInt     N;         
    PetscInt     N_elem_box;        
    PetscInt **box; 
    PetscScalar Lx, Ly, Lz; 
    PetscInt    Nx, Ny, Nz; 
    PetscScalar dx, dy, dz; 
    PetscScalar Xmin, Xmax; 
    PetscScalar Ymin, Ymax; 
    PetscScalar Zmin, Zmax; 
  } Boxes_Struct; 