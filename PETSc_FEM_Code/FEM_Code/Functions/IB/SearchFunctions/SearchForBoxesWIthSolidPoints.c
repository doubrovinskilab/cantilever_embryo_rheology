#ifndef SearchForTheBoxWithSolidPoints_C
#define SearchForTheBoxWithSolidPoints_C


PetscBool checkIfPointIn2DTriangle( Nodes_Struct *nodes, Elements_Struct *elements, Solid_Struct *solid, const PetscInt p, const PetscInt e)
{
  const PetscScalar xp = solid->x[p];
  const PetscScalar yp = solid->y[p];

  const PetscInt n0 = elements->inter[e][0];
  const PetscInt n1 = elements->inter[e][1];
  const PetscInt n2 = elements->inter[e][2];
  
  const PetscScalar x0 = nodes->x[n0];
  const PetscScalar y0 = nodes->y[n0];
  const PetscScalar x1 = nodes->x[n1];
  const PetscScalar y1 = nodes->y[n1];
  const PetscScalar x2 = nodes->x[n2];
  const PetscScalar y2 = nodes->y[n2];

  PetscScalar s = (-(x0 - x2)*(y0 - yp) + (x0 - xp)*(y0 - y2))/((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1));
  PetscScalar t = ( (x0 - x1)*(y0 - yp) - (x0 - xp)*(y0 - y1))/((x0 - x1)*(y0 - y2) - (x0 - x2)*(y0 - y1));

  return (((s>=0) && (s<=1) && (t>=0) && (t<=1) && ((s+t)<=1))? PETSC_TRUE: PETSC_FALSE);
}

PetscScalar Determinateof4x4(const PetscScalar a, const PetscScalar b, const PetscScalar c, const PetscScalar d,
  const PetscScalar e, const PetscScalar f, const PetscScalar g, const PetscScalar h,
  const PetscScalar i, const PetscScalar j, const PetscScalar k, const PetscScalar l,
  const PetscScalar m, const PetscScalar n, const PetscScalar o, const PetscScalar p)
{
  
  return a*(f*k*p - f*l*o - g*j*p + g*l*n + h*j*o - h*k*n) - b*(e*k*p - e*l*o - g*i*p + g*l*m + h*i*o - h*k*m) + c*(e*j*p - e*l*n - f*i*p + f*l*m + h*i*n - h*j*m) - d*(e*j*o - e*k*n - f*i*o + f*k*m + g*i*n - g*j*m);
}

PetscBool checkIfPointIn3DTetrahderon( Nodes_Struct *nodes, Elements_Struct *elements, Solid_Struct *solid, const PetscInt p, const PetscInt e)
{

  const PetscScalar xp = solid->x[p];
  const PetscScalar yp = solid->y[p];
  const PetscScalar zp = solid->z[p];

  const PetscInt n0 = elements->inter[e][0];
  const PetscInt n1 = elements->inter[e][1];
  const PetscInt n2 = elements->inter[e][2];
  const PetscInt n3 = elements->inter[e][3];

  const PetscScalar x0 = nodes->x[n0];
  const PetscScalar y0 = nodes->y[n0];
  const PetscScalar z0 = nodes->z[n0];
  const PetscScalar x1 = nodes->x[n1];
  const PetscScalar y1 = nodes->y[n1];
  const PetscScalar z1 = nodes->z[n1];
  const PetscScalar x2 = nodes->x[n2];
  const PetscScalar y2 = nodes->y[n2];
  const PetscScalar z2 = nodes->z[n2];
  const PetscScalar x3 = nodes->x[n3];
  const PetscScalar y3 = nodes->y[n3];
  const PetscScalar z3 = nodes->z[n3];

  const PetscScalar D0 = Determinateof4x4(x0, y0, z0, 1.0, x1, y1, z1, 1.0, x2, y2, z2, 1.0, x3, y3, z3, 1.0);
  const PetscScalar D1 = Determinateof4x4(xp, yp, zp, 1.0, x1, y1, z1, 1.0, x2, y2, z2, 1.0, x3, y3, z3, 1.0);
  const PetscScalar D2 = Determinateof4x4(x0, y0, z0, 1.0, xp, yp, zp, 1.0, x2, y2, z2, 1.0, x3, y3, z3, 1.0);
  const PetscScalar D3 = Determinateof4x4(x0, y0, z0, 1.0, x1, y1, z1, 1.0, xp, yp, zp, 1.0, x3, y3, z3, 1.0);
  const PetscScalar D4 = Determinateof4x4(x0, y0, z0, 1.0, x1, y1, z1, 1.0, x2, y2, z2, 1.0, xp, yp, zp, 1.0);

  return (((D0>=-EPS && D1>=-EPS && D2>=-EPS && D3>=-EPS && D4>=-EPS) || (D0<=EPS && D1<=EPS && D2<=EPS && D3<=EPS && D4<=EPS) )? PETSC_TRUE: PETSC_FALSE);
}

typedef PetscBool (*SearchForCheckingFunction)(Nodes_Struct *nodes, Elements_Struct *elements, Solid_Struct *solid,
  const PetscInt p, const PetscInt e); 

SearchForCheckingFunction DefineSearchForCheckingFunction(const PetscInt Dim)
{
  if(Dim == 2)
    return checkIfPointIn2DTriangle;
  else
    return checkIfPointIn3DTetrahderon;
}

PetscBool CheckIfPointInBox(Nodes_Struct *nodes, Elements_Struct *elements, 
  Boxes_Struct *boxes, Solid_Struct *solid, const PetscInt p, const PetscInt B, 
  const PetscMPIInt rank, SearchForCheckingFunction *CheckIfPointInElem)
{
  PetscBool ElemFound = PETSC_FALSE;
  
  for (PetscInt i = 0; i < boxes->box[B][0]; ++i) 
  {
    const PetscInt e = boxes->box[B][i+1]; 
    ElemFound = (*CheckIfPointInElem)(nodes, elements, solid, p, e);
    if ( ElemFound )  
    {
      solid->fluid_elem[p] = e, solid->N_total_pts_loc++; 
      solid->owner[p] = rank;
      break; 
    }
  }
  return ElemFound;
}

PetscBool CheckNeighbours(Nodes_Struct *nodes, Elements_Struct *elements,
  Boxes_Struct *boxes, Solid_Struct *solid, const PetscInt p, const PetscInt N_neig, 
  const PetscInt Dim, const PetscInt rank, SearchForCheckingFunction *CheckIfPointInElem)
{
  PetscInt i,j,k, B;
  PetscBool ElemFound = PETSC_FALSE;
  for(i = (-1*N_neig); i < (N_neig+1); i++)
  {
    for(j = (-1*N_neig); j < (N_neig+1); j++)
    {
      for(k = (-1*N_neig); k < (N_neig+1); k++)
      {
        B = BoxesCalcIndex(boxes, solid->x[p], solid->y[p], solid->z[p], i, j, k, Dim); 
        ElemFound = CheckIfPointInBox(nodes, elements, boxes, solid, p, B, rank, CheckIfPointInElem);
        if (ElemFound) break;
      }
      if (ElemFound) break;
    }
    if (ElemFound) break;
  }

  return ElemFound;
}


PetscBool CheckIfAllPointsHaveElements(Solid_Struct *solid)
{
  PetscBool AllPointsHaveElements = PETSC_TRUE;
  for (PetscInt p = 0; p < solid->N_total_pts; ++p) 
    if( solid->owner[p] == -1)
    { 
      AllPointsHaveElements = PETSC_FALSE;
      break;
    }

  return AllPointsHaveElements;
}

PetscErrorCode SearchForBoxWithPoints(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements, Boxes_Struct *boxes, Solid_Struct *solid)
{
  PetscLogDouble Clock1, Clock2;
  PetscTime(&Clock1);

  const PetscInt Dim  = system->Dim;
  const PetscInt rank = system->rank;
  PetscInt B;
  PetscBool ElemFound = PETSC_FALSE; 

  
  SearchForCheckingFunction CheckIfPointInElem  = DefineSearchForCheckingFunction(Dim);

  PetscPrintf(PETSC_COMM_WORLD,"  Searching for elements in boxes that hold the solid on all ranks\n");

  for(PetscInt j = 0; j < solid->N_total_pts; ++j) solid->owner[j] = -1;
  solid->N_total_pts_loc = 0; 

  for (PetscInt p = 0; p < solid->N_total_pts; ++p) 
  {
    ElemFound = PETSC_FALSE; 

    
    if (solid->fluid_elem[p] != -1) 
    { 
      ElemFound = CheckIfPointInElem(nodes, elements, solid, p, solid->fluid_elem[p]);
      if ( ElemFound ) 
      {
        solid->N_total_pts_loc++; 
        solid->owner[p] = system->rank;
        continue; 
      }
    }

    
    solid->fluid_elem[p] = -1;      

    B = BoxesCalcIndex(boxes, solid->x[p], solid->y[p], solid->z[p], 0, 0, 0, Dim); 
    ElemFound = CheckIfPointInBox(nodes, elements, boxes, solid, p, B, rank, &CheckIfPointInElem);
  }

  SolidBroadCastOwners(system, solid);

  if PetscNot(CheckIfAllPointsHaveElements(solid))
  {
    PetscPrintf(PETSC_COMM_WORLD,"    Checking neighbours on all ranks\n");

    
    for (PetscInt p = 0; p < solid->N_total_pts; ++p) 
    {
      ElemFound = PETSC_FALSE; 
      if( solid->owner[p] == -1) 
        ElemFound = CheckNeighbours(nodes, elements, boxes, solid, p, 1, Dim, rank, &CheckIfPointInElem);
    }

    SolidBroadCastOwners(system, solid);
    if PetscNot(CheckIfAllPointsHaveElements(solid))
    {
      
      PetscInt CFN = 0; 
      for (PetscInt p = 0; p < solid->N_total_pts; ++p) 
      {
        if( solid->owner[p] == -1)
        { 
          PetscErrorPrintf("Problem with point %d (%lf, %lf, %lf).\n", p, solid->x[p], solid->y[p], solid->z[p]); 
          if( BoxesCheckIfPointInDomain(boxes, solid->x[p], solid->y[p], solid->z[p], Dim) )
          { B = BoxesCalcIndex(boxes, solid->x[p], solid->y[p], solid->z[p], 0, 0, 0, Dim); 
            PetscErrorPrintf("No element was found for point %d which is in box %d with %d elements on rank %d\n", p, B, boxes->box[B][0], rank);   }
          else
          {PetscErrorPrintf("Point %d is not in the box domain of rank %d.\n", p, rank);    }
          SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error element was not found for point");     
        }
        if(solid->fluid_elem[p] != -1) CFN++;
      }

      if (CFN != solid->N_total_pts_loc)
      { PetscErrorPrintf("Number of points found do not match! Points with elem != -1 is %d while they should be %d\n", CFN, solid->N_total_pts_loc);
        SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Not all points are found"); }

    }
  }
  PetscBarrier(NULL); 
  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"    Finished searching in %lf seconds\n", (Clock2 - Clock1) );

  return 0;
}






#endif