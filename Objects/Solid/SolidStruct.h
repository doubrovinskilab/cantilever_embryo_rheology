typedef struct  
{
    SOLID_FORCE_LE_TYPE   Type;        
    PetscScalar           K;           
    PetscScalar           L0;          
    PetscScalar          *M;           
    PetscInt              M_size;      
} Forces_LE_Struct;

typedef struct  
 {
    SOLID_FORCE_SF_TYPE     Type;   
    PetscInt                  N;    
    PetscInt                 *p;    
    PetscScalar         A, B, C;    
    PetscScalar          Bd, Bv;    
    PetscScalar              ds;    
    PetscScalar              sf;    
} Forces_SF_Struct; 


typedef struct  
{
    PetscInt      N;                
    PetscInt     *p, *g;            
    PetscScalar  *x, *y, *z;          
    char          File[256];        
} Forces_BC_Speed_Struct;        

typedef struct  
{
    PetscInt      N;                
    PetscInt      *p, *g;           
    PetscScalar   *theta;           
    char          File[256];      
} Forces_BC_Speed_Ellipse_Struct;        

typedef struct 
{
    PetscScalar  *Fx, *Fy, *Fz; 
    PetscScalar  *x0, *y0, *z0; 
    Forces_LE_Struct    LE;      
    Forces_SF_Struct    SF;        
    Forces_BC_Speed_Struct          BC_FS; 
    Forces_BC_Speed_Ellipse_Struct  BC_FSE;
    PetscScalar   c;    
} Forces_Struct; 

typedef struct 
{
    SOLID_CONSTRAIN_TYPE   Type;    
    PetscInt                  N;    
    PetscInt                 *p;    
    PetscScalar         A, B, C;    
    PetscScalar          Bd, Bv;    
    PetscInt             axis_i;    
    PetscScalar          axis_v;    
    PetscScalar              ds;    
} Constrain_Struct; 


typedef struct 
{
    PetscInt     N_total_pts;    
    PetscInt     N_total_lns;    
    PetscInt     N_total_trs;    
    PetscInt     N_total_pts_loc; 
    PetscInt     Dim;             
    PetscScalar  *x,  *y,  *z;   
    PetscScalar  *Ux, *Uy, *Uz;  
    PetscScalar  *Sp;            
    PetscScalar  *Nx, *Ny, *Nz;  
    PetscInt     *fluid_elem;    
    PetscInt     *owner;         
    PetscInt    **lns;           
    PetscScalar  *ln_size;       
    PetscInt    **trs;           
    PetscScalar  *tr_size;           
    PetscInt    N_bound;            
    PetscInt    N_bound_pts;        
    PetscInt    N_bound_lns;        
    PetscInt    N_bound_trs;        
    PetscInt    *bound_pts;         
    PetscInt    *bound_lns;         
    PetscInt    *bound_trs;         
    PetscInt    *bound_group_pts;   
    PetscInt    *bound_group_lns;   
    PetscInt    *bound_group_trs;   
    PetscInt                    N_Gr;  
    SOLID_BC_TYPE           *Gr_type;  
    SOLID_ELEMENT_TYPE *Gr_elem_type;  
    PetscScalar           **Gr_value;  
    Forces_Struct         forces;
    Constrain_Struct      constrain;
} Solid_Struct; 