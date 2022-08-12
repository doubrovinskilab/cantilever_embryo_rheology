
// Function Type
typedef enum {FUNC_CONSTANT=0, FUNC_LINEAR=1, FUNC_EXP=2} FUNCTION_TYPE;
const char* FUNCTION_TYPE_CHAR[] = {"Constant", "Linear", "Exponential"};

// System Type
typedef enum {STOKES=0, IBM=1, DIFF=2, LINEAR_ELASTICITY=3} SYSTEM_TYPE;
const char* SYSTEM_TYPE_CHAR[] = {"STOKES", "IBM", "DIFFUSION", "LINEAR_ELASTICITY"};

// Fluid Solver Type
typedef enum { IMPLICIT=0, EXPLICIT=1, EXPLICIT_LUMPED=2, IMPLICIT_EXPLICIT_LUMPED=3} FLUID_SOLVER_TYPE;
const char* FLUID_SOLVER_TYPE_CHAR[] = {"Implicit", "Explicit", "Explicit + Lumped sum", "Implicit/Explicit + Lumped sum"};

// Fluid Solver Method
typedef enum { DIRECT=0, NEST=1, MONO=2, MONO_PRE=3} FLUID_SOLVER_METHOD;
const char* FLUID_SOLVER_METHOD_CHAR[] = {"Direct method", "Nested method", "Monolithic method", "Monolithic method with a Preconditioner"};

// Mesh Type
typedef enum {GMSH=0, VTK=1, FOAM=2, GMSH_OLD=3} MESH_TYPE;
const char* MESH_TYPE_CHAR[] = {"GMSH", "VTK", "FOAM", "GMSH_OLD"};

// Fluid Element Type
typedef enum { NO_ET=0, P0=1, P1=2, P2=3, P1b=4, P2b=5, Q1=6, Q2=7, Q2b=8} FLUID_ELEMENT_TYPE; 
const char* FLUID_ELEMENT_TYPE_CHAR[] = {"NO_ET", "P0", "P1", "P2", "P1b", "P2b", "Q1", "Q2", "Q2b"};

// Solid Element Type
typedef enum {NO_SET=0, POINT=1, LINE=2, TRIANGLE=3} SOLID_ELEMENT_TYPE; 
const char* SOLID_ELEMENT_TYPE_CHAR[] = {"NO_SET", "POINT", "LINE", "TRIANGLE"};

// Fluid Boundary Condition Type
typedef enum {INTERNAL=0, DIRICHLET=1, NEUMANN=2, ROBIN=3, SLIP=4, ZERO_NORMAL_VELOCITY=5, EXP_SLIP=6} FLUID_BC_TYPE;
const char* FLUID_BC_TYPE_CHAR[] = {"INTERNAL", "DIRICHLET", "NEUMANN", "ROBIN", "SLIP", "ZERO_NORMAL_VELOCITY", "EXP_SLIP"};

// Solid Boundary Condition Type
typedef enum {NO_SOLID_BC=0, FORCE_ADD=1, FORCE_INSERT=2, DIST_FIXED=3, DIST_SPEED=4, FORCE_SPEED=5, SYMMETRY=6 , 
  FORCE_SPEED_ELLIPSE=7, DIST_SPEED_ELLIPSE=8, DIST_SPEED_X_LINEAR_Y=9, AXIS_FIXED=10, FORCE_ADD_ELLIPSE=11} SOLID_BC_TYPE;
const char* SOLID_BC_TYPE_CHAR[] = {"NO_SOLID_BC", "FORCE_ADD", "FORCE_INSERT", "DIST_FIXED", "DIST_SPEED", 
  "FORCE_SPEED", "SYMMETRY", "FORCE_SPEED_ELLIPSE", "DIST_SPEED_ELLIPSE", "DIST_SPEED_X_LINEAR_Y", "AXIS_FIXED", "FORCE_ADD_ELLIPSE"};
/* Info about the SOLID BC:
  1- NO_SOLID_BC: No BC
  2- FORCE:       Input fixed boundary force at boundary points
  3- DIST_FIXED:  Fix the distance (do not move boundary points)
  4- DIST_SPEED:  Move the points by a prescribed boundary speed
  5- FORCE_SPEED: Create fake points and have them be moved at a prescribed boundary speed. Attach 
                  to each of these fake points a boundary point with an elastic connection. and have it be 
                  dragged behind theses fake point by a large elastic spring. 
  6- SYMMETRY:    Tries to enforce symmetry using one of the axis as a symmetry plane
  7- FORCE_SPEED_ELLIPSE: Create fake point along an ellipse and attach them to boundary points with an elastic spring.
                          Next step move the fake point and transfer the force to the boundary points 
  8- DIST_SPEED_ELLIPSE:  Move boundary points along an ellipse in a constant speed
  9- DIST_SPEED_X_LINEAR_Y: Move boundary points in the x direction with a velocity that varies linearly with y
  10- AXIS_FIXED:           Fix the position of one axis of the boundary point
  11- FORCE_ADD_ELLIPSE: Add a force that is tanget to the surface of an ellipse
*/

// Fluid Viscosity Models
typedef enum {ONE_CONSTANT=0, TWO_CONSTANTS_ONE_ELLIPSOID=1, TWO_CONSTANTS_TWO_ELLIPSOIDS=2} FLUID_VISCOSITY_METHOD;
const char* FLUID_VISCOSITY_METHOD_CHAR[] = {"ONE_CONSTANT", "TWO_CONSTANTS_ONE_ELLIPSOID", "TWO_CONSTANTS_TWO_ELLIPSOIDS"};
/* Info about the FLUID_VISCOSITY_METHOD:
  1- SINGLE_CONSTANT:                     Only one single constant visocisuty for all the points
  2- TWO_CONSTANTS_ELLIPSOID:             Location of a point with respect to an ellipsoid will determin which of the two viscosities it will take
  3- TWO_CONSTANTS_TWO_ELLIPSOIDS:        Location of a point with respect to two ellipsoids will determin which of the two viscosities it will take.
                                          An ellipsoid (A, Bd, C) for y>0 and another ellipsoid (A, Bv, C) for y<=0*/

// Constrain: 
typedef enum {CONSTRAIN_NULL=0, CONSTRAIN_ELLIPSOID=1, CONSTRAIN_TWO_ELLIPSOIDS=2, CONSTRAIN_AXIS=3} SOLID_CONSTRAIN_TYPE;
const char* SOLID_CONSTRAIN_CHAR[] = {"NO_CONSTRAINT", "CONSTRAIN_ELLIPSOID", "CONSTRAIN_TWO_ELLIPSOIDS", "CONSTRAIN_AXIS"};
/* Info about the FLUID_VISCOSITY_METHOD:
  0- CONSTRAIN_NULL:                      No Geometric constrain on the model
  1- CONSTRAIN_ELLIPSOID:                 Constrain the point on one ellipsoid
  2- CONSTRAIN_TWO_ELLIPSOIDS:            Constrain the point on two ellipsoids [An ellipsoid (A, Bd, C) for y>0 and another ellipsoid (A, Bv, C) for y<=0]
  3- CONSTRAIN_AXIS:                      Constrain on one axis
  */

// L-Node Dynamics:
typedef enum {LNODE_DYNAMICS_NULL=0, LNODE_DYNAMICS_ALWAYS=1, LNODE_DYNAMICS_COND=2} SOLID_LNODE_DYNAMICS_TYPE;
const char* SOLID_LNODE_DYNAMICS_CHAR[] = {"NO_L-NODE_DYNAMICS", "L-LNODE_DYNAMICS_ALWAYS", "L-NODE_DYNAMICS_COND"};
/* Info about the FLUID_VISCOSITY_METHOD:
  0- LNODEDYNAMICS_NULL:     No L-Node Dynamics
  1- LNODEDYNAMICS_ALAWAYS:  L-Node Dynamics are always one (i.e. L0_new = dt/tau*(L-L0_old) + L0_old)  [NOT DEFINED]
  2- LNODEDYNAMICS_COND:     L-Node Dynamics are edge dependent and they are activated if L > L0*lambda [NOT DEFINED]
  */

// Noise Models
typedef enum {NOISE_NULL=0, NOISE_CONST=1} SOLID_NOISE_TYPE;
const char* SOLID_NOISE_CHAR[] = {"NO_NOISE_ADDED", "NOISE_CONST"};
/* Info about the SOLID_FORCE_ST_TYPE:
  0- NOISE_NULL:      No noise
  1- NOISE_CONST:     Const noise*/

// Forces: Linear Elasticity Models
typedef enum {FORCE_LE_NULL=0, FORCE_LE_FIXED_DISP=1, FORCE_LE_INITIAL_DISP=2, FORCE_LE_INITIAL_DISP_MULTI=3, FORCE_LE_INITIAL_DISP_COND=4} SOLID_FORCE_LE_TYPE;
const char* SOLID_FORCE_LE_CHAR[] = {"NO_LINEAR_ELASTIC_FORCE", "LINEAR_ELASTIC_FORCE_FIXED_LENGTH", "LINEAR_ELASTIC_FORCE_INITIAL_LENGTH", "LINEAR_ELASTIC_FORCE_MULTIPLE_INITIAL_LENGTH", 
  "LINEAR_ELASTIC_FORCE_INITIAL_LENGTH_WITH_CONDITION"};
/* Info about the Solid Internal Forces Type:
  0- FORCE_LE_NULL:                     No Forces between points
  1- FORCE_LE_FIXED_DISP:               Linear Elastic Force = K*(D-D0) where D0 is from input saved in system
  2- FORCE_LE_INITIAL_DISP:             Linear Elastic Force = K*(D-D0) where D0 is calculated at t=0 from initital configutation
  3- FORCE_LE_INITIAL_DISP_MULTI:       Linear Elastic Force = K*(D-D0) where D0 is calculated at t=0 from initital configutation &
                                          K varies at each spring according to a dat file 
  4- FORCE_LE_INITIAL_DISP_COND:        Linear Elastic Force = K*(D-D0) if D>D0 else F = 0, where D0 is calculated at t=0 from initital configutation [NOT  DEFINED]
                                  */

// Forces: Surface Tension Models
typedef enum {FORCE_ST_NULL=0, FORCE_ST_CONST=1, FORCE_ST_MULTI=2, FORCE_ST_CONST_PERC=3, FORCE_ST_MULTI_PERC=4} SOLID_FORCE_ST_TYPE;
const char* SOLID_FORCE_ST_CHAR[] = {"NO_SURFACE_TENSION_FORCE", "SURFACE_TENSION_CONST", "SURFACE_TENSION_MULT", "SURFACE_TENSION_PERCENTAGE", "SURFACE_TENSION_MULTI_PERCENTAGE"};
/* Info about the SOLID_FORCE_ST_TYPE:
  0- FORCE_ST_NULL:         No Surface Tension
  1- FORCE_ST_CONST:        Const Surface Tension [NOT DEFINED]
  2- FORCE_ST_MULTI:        Multiple Surface Tension for each edge  [NOT DEFINED]
  3- FORCE_ST_CONST_PERC:   Const Surface Tension percetange from K [NOT DEFINED]
  4- FORCE_ST_MULTI_PERC:   Multiple Surface Tension percetange from K for each edge  [NOT DEFINED]
  */

// Forces: Volume Conservation Models
typedef enum {FORCE_VC_NULL=0, FORCE_VC_CONST=1} SOLID_FORCE_VC_TYPE;
const char* SOLID_FORCE_VC_CHAR[] = {"NO_VOLUME_CONSERVATION_FORCE", "VOLUME_CONSERVATION_CONST"};
/* Info about the SOLID_FORCE_VC_TYPE:
  0- FORCE_VC_NULL:      No Volume Conservation
  1- FORCE_VC_CONST:     Const Volume Conservation [NOT DEFINED] 
*/

// Forces: Soft Force Models
typedef enum {FORCE_SF_NULL=0, FORCE_SF_ELLIPSOID=1, FORCE_SF_TWO_ELLIPSOIDS=2, FORCE_SF_ELLIPSOID_REMOVE_NORMAL_FORCE=3} SOLID_FORCE_SF_TYPE;
const char* SOLID_FORCE_SF_CHAR[] = {"NO_SOFT_FORCE", "FORCE_SF_ELLIPSOID", "FORCE_SF_TWO_ELLIPSOIDS" , "FORCE_SF_ELLIPSOID_REMOVE_NORMAL_FORCE"};
/* Info about the SOLID_FORCE_SF_TYPE:
  0- FORCE_SF_NULL:                           
  1- FORCE_SF_ELLIPSOID:                      Add a soft force to constraint the point on one ellipsoid
  2- FORCE_SF_TWO_ELLIPSOIDS:                 Remove normal force components of the chosen ellipsoid
  3- FORCE_SF_ELLIPSOID_REMOVE_NORMAL_FORCE:  Add a soft force on two ellipsoids to remove normal component of force
  */


// Forces: Barrier Force Models
typedef enum {FORCE_BF_NULL=0, FORCE_BF_ALL_POINTS=1, FORCE_BF_BOUNDARY_POINTS=2} SOLID_FORCE_BF_TYPE;
const char* SOLID_FORCE_BF_CHAR[] = {"NO_BARRIER_FORCE", "BARRIER_FORCE_ALL_POINTS", "BARRIER_FORCE_BOUNDARY_POINTS"};
/* Info about the SOLID_FORCE_BF_TYPE:
  0- FORCE_BF_NULL:               No Barrier Force
  1- FORCE_BF_ALL_POINTS:         Barrier Force on all points  [NOT DEFINED] 
  2- FORCE_BF_BOUNDARY_POINTS:    Barrier Force on boundary points [NOT DEFINED] 
  */