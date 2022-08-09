#ifndef Boxes_Functions_C
#define Boxes_Functions_C

void BoxesInitialize(Boxes_Struct *boxes)
{
  boxes->N = -1;
  boxes->N_elem_box = -1;
  boxes->Lx = -1;
  boxes->Ly = -1;
  boxes->Lz = -1;
  boxes->Nx = -1;
  boxes->Ny = -1;
  boxes->Nz = -1;
  boxes->dx = -1;
  boxes->dy = -1;
  boxes->dz = -1;
  boxes->Xmin =  1e10;  
  boxes->Xmax = -1e10;  
  boxes->Ymin =  1e10;  
  boxes->Ymax = -1e10;  
  boxes->Zmin =  1e10;  
  boxes->Zmax = -1e10;  
}

PetscInt BoxesCalcIndex(Boxes_Struct *boxes, const PetscScalar Px, const PetscScalar Py, const PetscScalar Pz, 
  const PetscInt Ax, const PetscInt Ay, const PetscInt Az, const PetscInt Dim) 
{
  PetscInt ix = (PetscInt) ((Px - boxes->Xmin)/boxes->dx) + Ax; 
  PetscInt iy = (PetscInt) ((Py - boxes->Ymin)/boxes->dy) + Ay; 
  PetscInt iz = (PetscInt) ((Pz - boxes->Zmin)/boxes->dz) + Az; 
  ix = PetscMax(PetscMin(ix,boxes->Nx-1),0);
  iy = PetscMax(PetscMin(iy,boxes->Ny-1),0);
  iz = PetscMax(PetscMin(iz,boxes->Nz-1),0);
  return ix + iy*boxes->Nx + iz*boxes->Nx*boxes->Ny*(Dim-2); 
}

PetscBool BoxesCheckIfPointInDomain(Boxes_Struct *boxes, const PetscScalar Px, const PetscScalar Py, 
  const PetscScalar Pz, const PetscInt Dim) 
{
  PetscBool NodeInDomain = PETSC_TRUE; 
  if      (  Px < boxes->Xmin || Px > boxes->Xmax ) 
    NodeInDomain = PETSC_FALSE;
  else if (  Py < boxes->Ymin || Py > boxes->Ymax ) 
    NodeInDomain = PETSC_FALSE;
  else if ( (Pz < boxes->Zmin || Pz > boxes->Zmax ) && (Dim == 3) )
    NodeInDomain = PETSC_FALSE;
  return NodeInDomain;
}

PetscBool BoxesCheckIfElementInBox(Boxes_Struct *boxes, const PetscInt B, const PetscInt e)
{
  PetscBool ElementFound = PETSC_FALSE; 
  for (PetscInt k = 0; k < boxes->box[B][0]; k++) 
    if ( boxes->box[B][k+1] == e)
    { ElementFound = PETSC_TRUE;
      break;}
  return ElementFound;
}

PetscErrorCode BoxesMallocArrays(Boxes_Struct *boxes)
{
  PetscInt i,j;
  PetscPrintf(PETSC_COMM_WORLD," Creating the box array\n");
  PetscMalloc1(boxes->N, &boxes->box); 
  for(i = 0; i < boxes->N; ++i)
  { PetscMalloc1(boxes->N_elem_box + 1, &boxes->box[i]);
    boxes->box[i][0] = 0;
    for(j = 0; j < boxes->N_elem_box; ++j)  boxes->box[i][j+1] = -1;}
  PetscBarrier(NULL); 
  return 0;
}

PetscErrorCode BoxesFree(Boxes_Struct *boxes)
{
  for(PetscInt i = 0; i <boxes->N; ++i)  PetscFree(boxes->box[i]);
  PetscFree(boxes->box);
  return 0;
}


PetscErrorCode BoxesFindDomainParameters(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements, Boxes_Struct *boxes)
{
  PetscErrorCode ierr;
  for (PetscInt i = 0; i < nodes->N_U_total_loc; i++)
  {
    const PetscScalar x = nodes->x[i];
    const PetscScalar y = nodes->y[i];
    const PetscScalar z = nodes->z[i];

    if ( x > boxes->Xmax) boxes->Xmax = x;
    if ( y > boxes->Ymax) boxes->Ymax = y;
    if ( z > boxes->Zmax) boxes->Zmax = z;

    if ( x < boxes->Xmin) boxes->Xmin = x;
    if ( y < boxes->Ymin) boxes->Ymin = y;
    if ( z < boxes->Zmin) boxes->Zmin = z;
  }

  boxes->Lx = (boxes->Xmax - boxes->Xmin); 
  boxes->Ly = (boxes->Ymax - boxes->Ymin); 
  boxes->Lz = (boxes->Zmax - boxes->Zmin); 

  
  boxes->Xmin = boxes->Xmin - 0.025*boxes->Lx;
  boxes->Ymin = boxes->Ymin - 0.025*boxes->Ly;
  boxes->Zmin = boxes->Zmin - 0.025*boxes->Lz;
  boxes->Xmax = boxes->Xmax + 0.025*boxes->Lx;
  boxes->Ymax = boxes->Ymax + 0.025*boxes->Ly;
  boxes->Zmax = boxes->Zmax + 0.025*boxes->Lz;

  boxes->Lx = (boxes->Xmax - boxes->Xmin); 
  boxes->Ly = (boxes->Ymax - boxes->Ymin); 
  boxes->Lz = (boxes->Zmax - boxes->Zmin);

  const PetscInt n0 = elements->inter[0][0];
  const PetscInt n1 = elements->inter[0][1];
  const PetscScalar Dx = nodes->x[n0]-nodes->x[n1];
  const PetscScalar Dy = nodes->y[n0]-nodes->y[n1];
  const PetscScalar Dz = nodes->z[n0]-nodes->z[n1];

  const PetscScalar D = PetscSqrtScalar( Dx*Dx + Dy*Dy + Dz*Dz);  
  boxes->Nx = PetscMax( (PetscInt) (boxes->Lx/D - 1), 1);
  boxes->Ny = PetscMax( (PetscInt) (boxes->Ly/D - 1), 1);
  boxes->Nz = PetscMax( (PetscInt) (boxes->Lz/D - 1), 1);
  if (system->Dim == 2) boxes->Nz = 1;
  boxes->N = boxes->Nx*boxes->Ny*boxes->Nz;
  boxes->dx = boxes->Lx / boxes->Nx;
  boxes->dy = boxes->Ly / boxes->Ny;
  boxes->dz = boxes->Lz / boxes->Nz;

  ierr = PetscBarrier(NULL); CHKERRQ(ierr);
  return 0;
}

PetscErrorCode BoxesCountMaxElements(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements, Boxes_Struct *boxes)
{  
  const PetscInt N_iE = elements->N_inter_loc;
  const PetscInt P_iS = elements->P_S_inter;

  PetscScalar Ex, Ey, Ez;
  PetscInt B;
  PetscInt *Count;
  PetscCalloc1(boxes->N, &Count);
  boxes->N_elem_box = 0;
  for (PetscInt e = 0; e < N_iE; e++) 
  {
    CalcElementCentroid(nodes->x, nodes->y, nodes->z, elements->inter[e], elements->P_S_inter, &Ex, &Ey, &Ez);
    B = BoxesCalcIndex(boxes, Ex, Ey, Ez, 0, 0, 0, system->Dim);
    Count[B]++; 
    if (Count[B] > boxes->N_elem_box) boxes->N_elem_box = Count[B];
  }
  for (PetscInt e = 0; e < N_iE; e++) 
    for (PetscInt m = 0; m < P_iS; m++) 
    {
        const PetscInt    n = elements->inter[e][m];
        const PetscScalar X = nodes->x[n];
        const PetscScalar Y = nodes->y[n];
        const PetscScalar Z = nodes->z[n];
        if (BoxesCheckIfPointInDomain(boxes, X, Y, Z, system->Dim) )
        {
          B = BoxesCalcIndex(boxes, X, Y, Z, 0, 0, 0, system->Dim);
          Count[B]++; 
          if (Count[B] > boxes->N_elem_box) boxes->N_elem_box = Count[B];
        }
    }
  PetscFree(Count);
  PetscBarrier(NULL); 

  return 0;
}


PetscErrorCode BoxesPlaceElements(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements, Boxes_Struct *boxes)
{
  PetscScalar Ex, Ey, Ez;
  PetscInt B;
  const PetscInt  N_iE = elements->N_inter_loc;
  const PetscInt  P_iS = elements->P_S_inter;
  for (PetscInt e = 0; e < N_iE; e++) 
  {
    CalcElementCentroid(nodes->x, nodes->y, nodes->z, elements->inter[e], elements->P_S_inter, &Ex, &Ey, &Ez);
    B = BoxesCalcIndex(boxes, Ex, Ey, Ez, 0, 0, 0, system->Dim);
    boxes->box[B][0]++; 
    boxes->box[B][boxes->box[B][0]] = e;
  }
  PetscBarrier(NULL); 
  for (PetscInt e = 0; e < N_iE; e++) 
    for (PetscInt m = 0; m < P_iS; m++) 
    {
      { 
        const PetscInt    n = elements->inter[e][m];
        const PetscScalar X = nodes->x[n];
        const PetscScalar Y = nodes->y[n];
        const PetscScalar Z = nodes->z[n];
        if (BoxesCheckIfPointInDomain(boxes, X, Y, Z, system->Dim) )
        {
          B = BoxesCalcIndex(boxes, X, Y, Z, 0, 0, 0, system->Dim);
          if ( PetscNot(BoxesCheckIfElementInBox(boxes, B, e)) )
          { boxes->box[B][0]++; 
            if ( boxes->box[B][0] > boxes->N_elem_box) {  
              SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG, "Number of elements in box exceeded the max.");  }  
            boxes->box[B][boxes->box[B][0]] = e; 
          }
        }
      }
    }
  PetscBarrier(NULL); 

  return 0;
}


PetscErrorCode BoxesCreateDomain(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements, Boxes_Struct *boxes)
{
  BoxesFindDomainParameters(system, nodes, elements, boxes);
  BoxesCountMaxElements(system, nodes, elements, boxes);
  BoxesMallocArrays(boxes);
  BoxesCountMaxElements(system, nodes, elements, boxes);
  BoxesPlaceElements(system, nodes, elements, boxes);
  return 0;
}

#endif 