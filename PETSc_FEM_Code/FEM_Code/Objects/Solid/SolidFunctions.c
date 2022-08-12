#ifndef Solid_Functions_C
#define Solid_Functions_C

#include "ConstrainFunctions.c"

void SolidInitialize(Solid_Struct *solid)
{
  solid->N_total_pts = -1;
  solid->N_total_lns = -1;
  solid->N_total_trs = -1;
  solid->N_total_pts_loc = -1;
  solid->N_bound = -1;
  solid->N_bound_pts = -1;
  solid->N_bound_lns = -1;
  solid->N_bound_trs = -1;
  solid->N_Gr = -1;
}

void SolidInitializeProperties(System_Struct *system, Solid_Struct *solid)
{ 
  ConstrainCreate(system, solid);
}

void SolidAddProperties(System_Struct *system, Solid_Struct *solid)
{
}


PetscErrorCode SolidMallocPointArrays(Solid_Struct *solid) 
{
  if(solid->N_total_pts <= 0)  {
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error solid->N_total_pts");      }
  
  PetscMalloc3(solid->N_total_pts, &solid->x,  solid->N_total_pts, &solid->y,  solid->N_total_pts, &solid->z);
  PetscMalloc3(solid->N_total_pts, &solid->Ux, solid->N_total_pts, &solid->Uy, solid->N_total_pts, &solid->Uz);
  PetscMalloc3(solid->N_total_pts, &solid->Nx, solid->N_total_pts, &solid->Ny, solid->N_total_pts, &solid->Nz);
  PetscMalloc3(solid->N_total_pts, &solid->Sp, solid->N_total_pts, &solid->fluid_elem, solid->N_total_pts, &solid->owner);
  for(PetscInt i=0; i<solid->N_total_pts; i++)
    solid->Sp[i]=-1, solid->fluid_elem[i] = -1, solid->owner[i]=-1;

  return 0; 
}

PetscErrorCode SolidMallocCellArrays(Solid_Struct *solid) 
{
  if(solid->N_total_lns < 0 || solid->N_total_trs < 0)  {
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error solid->N_total_lns & solid->N_total_trs");      }
  
  PetscMalloc1(solid->N_total_lns, &solid->lns);
  for (PetscInt i=0; i<solid->N_total_lns; i++) 
    PetscMalloc1(2, &solid->lns[i]), solid->lns[i][0] = -1, solid->lns[i][1] = -1;

  PetscMalloc1(solid->N_total_trs, &solid->trs);
  for (PetscInt i=0; i<solid->N_total_trs; i++) 
    PetscMalloc1(3, &solid->trs[i]), solid->trs[i][0] = -1, solid->trs[i][1] = -1, solid->trs[i][2] = -1;
  return 0; 
}

PetscErrorCode SolidMallocBoundaryCellArrays(Solid_Struct *solid) 
{
  if(solid->N_bound_pts<0 || solid->N_bound_pts<0 || solid->N_bound_trs<0)  {
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error solid->N_bound_pts or solid->N_bound_pts or solid->N_bound_trs");      }
  
  PetscMalloc2(solid->N_bound_pts, &solid->bound_pts, solid->N_bound_pts, &solid->bound_group_pts);
  PetscMalloc2(solid->N_bound_lns, &solid->bound_lns, solid->N_bound_lns, &solid->bound_group_lns);
  PetscMalloc2(solid->N_bound_trs, &solid->bound_trs, solid->N_bound_trs, &solid->bound_group_trs);

  return 0; 
}

PetscErrorCode SolidMallocBoundaryGroupArrays(Solid_Struct *solid) 
{
  if(solid->N_Gr<0 )  {
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error solid->N_Gr ");      }
  
  PetscMalloc2(solid->N_Gr, &solid->Gr_type, solid->N_Gr, &solid->Gr_elem_type);
  PetscMalloc1(solid->N_Gr, &solid->Gr_value);
  for (PetscInt i=0; i<solid->N_Gr;i++) {
    solid->Gr_elem_type[i]=NO_SET;
    PetscMalloc1(10, &solid->Gr_value[i]);  }

  return 0; 
}


PetscErrorCode SolidFree(Solid_Struct *solid) 
{
  PetscPrintf(PETSC_COMM_WORLD,"  Solid memory is being freed\n");

  PetscFree3(solid->x,  solid->y,  solid->z);
  PetscFree3(solid->Ux, solid->Uy, solid->Uz);
  PetscFree3(solid->Nx, solid->Ny, solid->Nz);
  PetscFree3(solid->Sp, solid->fluid_elem, solid->owner);

  for (PetscInt i=0; i<solid->N_total_lns;i++)  PetscFree(solid->lns[i]);
  PetscFree(solid->lns);

  for (PetscInt i=0; i<solid->N_total_trs;i++)  PetscFree(solid->trs[i]);
  PetscFree(solid->trs);

  PetscFree2(solid->bound_pts,  solid->bound_group_pts);
  PetscFree2(solid->bound_lns,  solid->bound_group_lns);
  PetscFree2(solid->bound_trs,  solid->bound_group_trs);

  PetscFree2(solid->Gr_type, solid->Gr_elem_type);
  for (PetscInt i=0; i<solid->N_Gr;i++) PetscFree(solid->Gr_value[i]);
  PetscFree(solid->Gr_value);

  ConstrainFree(solid);
  
  return 0; 
}

PetscErrorCode SolidPrintInfo(System_Struct *system, Solid_Struct *solid, PetscBool boolSolid) 
{
  if (boolSolid)
  {

    PetscPrintf(PETSC_COMM_WORLD,"  Total number of points:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d\n", system->rank, solid->N_total_pts);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"    x:\n");
    PetscScalarView(solid->N_total_pts, solid->x, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"    y:\n");
    PetscScalarView(solid->N_total_pts, solid->y, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"    z:\n");
    PetscScalarView(solid->N_total_pts, solid->z, stdout_viewer); 

    PetscPrintf(PETSC_COMM_WORLD,"  Total number of groups:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d\n", system->rank, solid->N_Gr);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

    PetscPrintf(PETSC_COMM_WORLD,"  Total number of N_total_lns:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d\n", system->rank, solid->N_total_lns);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"    lns:\n");
    for(PetscInt i=0; i<solid->N_total_lns;i++) PetscIntView(2, solid->lns[i], stdout_viewer);

    PetscPrintf(PETSC_COMM_WORLD,"  Total number of N_total_trs:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d\n", system->rank, solid->N_total_trs);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"    trs:\n");
    for(PetscInt i=0; i<solid->N_total_trs;i++) PetscIntView(3, solid->trs[i], stdout_viewer);

    PetscPrintf(PETSC_COMM_WORLD,"  Total number of N_bound:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d\n", system->rank, solid->N_bound);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);

    PetscPrintf(PETSC_COMM_WORLD,"  Total number of N_bound_pts:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d\n", system->rank, solid->N_bound_pts);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"    bound_pts:\n");
    PetscIntView(solid->N_bound_pts, solid->bound_pts, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"    bound_group_pts:\n");
    PetscIntView(solid->N_bound_pts, solid->bound_group_pts, stdout_viewer);

    PetscPrintf(PETSC_COMM_WORLD,"  Total number of N_bound_lns:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d\n", system->rank, solid->N_bound_lns);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"    bound_lns:\n");
    PetscIntView(solid->N_bound_lns, solid->bound_lns, stdout_viewer);
    PetscPrintf(PETSC_COMM_WORLD,"    bound_group_lns:\n");
    PetscIntView(solid->N_bound_lns, solid->bound_group_lns, stdout_viewer);

    PetscPrintf(PETSC_COMM_WORLD,"  Total number of N_bound_trs:\n");
    PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d:%d\n", system->rank, solid->N_bound_trs);
    PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT);
    PetscPrintf(PETSC_COMM_WORLD,"    bound_trs:\n");
    PetscIntView(solid->N_bound_trs, solid->bound_trs, stdout_viewer); 
    PetscPrintf(PETSC_COMM_WORLD,"    bound_group_trs:\n");
    PetscIntView(solid->N_bound_trs, solid->bound_group_trs, stdout_viewer); 
  }
  return 0;
}


PetscErrorCode SolidBroadCastOwners(System_Struct *system, Solid_Struct *solid)
{
  PetscErrorCode ierr;
  MPI_Request Request;
  MPI_Status  Status;
  PetscInt    j;
  PetscInt  size_s, *small_owner;

  PetscBarrier(NULL);
  if (system->size != 1) 
  {
    if (system->rank != 0) 
    { 

      size_s = 0;
      for(j= 0; j < solid->N_total_pts; ++j) 
        if (solid->owner[j] != -1) size_s++;

      PetscMalloc1(size_s, &small_owner);
      PetscInt i=0;
      for(j= 0; j < solid->N_total_pts; ++j) 
        if (solid->owner[j] != -1){ small_owner[i] = j; i++;}

      
      ierr = MPI_Isend(&size_s, 1, MPIU_INT, 0, system->rank+10000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

      ierr = MPI_Isend(small_owner, size_s, MPIU_INT, 0, system->rank+20000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

      PetscFree(small_owner);
    } 
    else
    {
      for(PetscInt r = 1; r<system->size; ++r) 
      {
        size_s = -1;
        ierr = MPI_Irecv(&size_s, 1, MPIU_INT,  r, r+10000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

        PetscMalloc1(size_s, &small_owner);
        ierr = MPI_Irecv( small_owner, size_s, MPIU_INT,  r, r+20000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

        for(j= 0; j < size_s; ++j) 
            solid->owner[small_owner[j]] = r;

        PetscFree(small_owner);
      }
    }
    PetscBarrier(NULL);
  
    
    ierr = MPI_Bcast(solid->owner, solid->N_total_pts, MPIU_INT, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
  }
  return 0;
}


PetscErrorCode SolidCalculateLineSize(Solid_Struct *solid)
{
  PetscMalloc1(solid->N_total_lns, &solid->ln_size);
  for (PetscInt j=0; j<solid->N_total_lns; j++){
    solid->ln_size[j] = CalcDistance(solid->x, solid->y, solid->z, solid->lns[j][0], solid->lns[j][1]);
  }
  return 0;
}


PetscErrorCode SolidCalculateTriangleSize(Solid_Struct *solid)
{
  PetscMalloc1(solid->N_total_trs, &solid->tr_size);
  for (PetscInt j=0; j<solid->N_total_trs; j++){
    solid->tr_size[j] = CalcAreaOfTriangle(solid->x, solid->y, solid->z, solid->trs[j][0], solid->lns[j][1], solid->trs[j][2]);
  }
  return 0;
}


PetscErrorCode SolidConvertTrs2Lns(System_Struct *system, Solid_Struct *solid)
{ 

  if (system->SDMesh_ConvTrs2Lns)
  {
    if (solid->N_total_trs != 0)
    {
      PetscInt N_lns=0;             
      PetscBool l_added=PETSC_FALSE;        
      PetscInt  N_lns_tmp = solid->N_total_trs*3+solid->N_total_lns; 
      PetscInt **lns_tmp;
      PetscMalloc1(N_lns_tmp, &lns_tmp);

      for(PetscInt j = 0; j < N_lns_tmp; ++j) PetscMalloc1(2, &lns_tmp[j]);
      
      for(PetscInt j = 0; j < solid->N_total_lns; ++j) 
        lns_tmp[j][0] = solid->lns[j][0], lns_tmp[j][1] = solid->lns[j][1];
      N_lns = solid->N_total_lns;

      
      PetscInt i, l;
      for (PetscInt t=0; t<solid->N_total_trs; t++)
        for(i=0; i <3; i++)
        {
          const PetscInt p0 = solid->trs[t][i];         
          const PetscInt p1 = solid->trs[t][(i+1)%3];   

          
          
          l_added=PETSC_FALSE;               
          for(l=0; l<N_lns; l++)
          { 
            
            if(lns_tmp[l][0] == p0 && lns_tmp[l][1] == p1){
                l_added = PETSC_TRUE;
                break;}

            if(lns_tmp[l][0] == p1 && lns_tmp[l][1] == p0){
                l_added = PETSC_TRUE;
                break;}
          }

          if (PetscNot(l_added))  
          {                 
            lns_tmp[N_lns][0] = p0;
            lns_tmp[N_lns][1] = p1;
            N_lns++;
          }    
        }

      
      for (PetscInt j=0; j<solid->N_total_lns;j++)  PetscFree(solid->lns[j]);
      PetscFree(solid->lns); 

      
      solid->N_total_lns = N_lns;
      PetscMalloc1(solid->N_total_lns, &solid->lns);
      for (PetscInt j=0; j<solid->N_total_lns;j++) {
        PetscMalloc1(2, &solid->lns[j]); 
        solid->lns[j][0] = lns_tmp[j][0]; 
        solid->lns[j][1] = lns_tmp[j][1];}

      
      for(PetscInt j = 0; j < N_lns_tmp; ++j) PetscFree(lns_tmp[j]);
      PetscFree(lns_tmp);
    }
  }
  return 0;
}

PetscErrorCode SolidCalculateWeightedCenter(Solid_Struct *solid, PetscScalar *Gx, 
  PetscScalar *Gy, PetscScalar *Gz)
{
  PetscScalar Total_S = 0;
  *Gx = 0, *Gy = 0, *Gz=0;
  for(PetscInt p = 0; p < solid->N_total_pts; ++p)
  {
    Total_S += solid->Sp[p];
    *Gx += solid->x[p]*solid->Sp[p];
    *Gy += solid->y[p]*solid->Sp[p];
    *Gz += solid->z[p]*solid->Sp[p];
  }
  *Gx = *Gx/Total_S;
  *Gy = *Gy/Total_S;
  *Gz = *Gz/Total_S;

  return 0;
}

PetscErrorCode SolidCalculateGeometricCenter(Solid_Struct *solid, PetscScalar *Gx, 
  PetscScalar *Gy, PetscScalar *Gz)
{
  *Gx = 0, *Gy = 0, *Gz=0;
  for(PetscInt p = 0; p < solid->N_total_pts; ++p)
  {
    *Gx += solid->x[p]*solid->Sp[p];
    *Gy += solid->y[p]*solid->Sp[p];
    *Gz += solid->z[p]*solid->Sp[p];
  }
  *Gx = *Gx/solid->N_total_pts;
  *Gy = *Gy/solid->N_total_pts;
  *Gz = *Gz/solid->N_total_pts;

  return 0;
}

PetscErrorCode SolidCalculateSizeAndNormalAtEachPoint(System_Struct *system, Solid_Struct *solid)
{
  for (PetscInt p = 0; p < solid->N_total_pts; ++p)
    solid->Sp[p] = 0, solid->Nx[p] = 0, solid->Ny[p] = 0, solid->Nz[p]=0;

  if (solid->Dim == 0) 
  {
    return 0;
  }
  else if (system->Dim == 2 && (solid->Dim == 2 || solid->Dim == 1)) 
  {

    PetscScalar *Tx, *Ty, *Tz;    
    PetscScalar  Vx,  Vy,  Vz;
    PetscScalar  Gx,  Gy,  Gz;    
    PetscScalar Vgpx, Vgpy, Vgpz; 

    
    PetscMalloc3(solid->N_total_pts, &Tx, solid->N_total_pts, &Ty, solid->N_total_pts, &Tz); 
    for(PetscInt p = 0; p < solid->N_total_pts; ++p) 
      Tx[p] = 0, Ty[p] = 0, Tz[p] = 0;

    
    for(PetscInt l = 0; l < solid->N_total_lns; ++l)    
    {
      const PetscInt p0 = solid->lns[l][0];  
      const PetscInt p1 = solid->lns[l][1];  
      CalcVector( solid->x, solid->y, solid->z, p0, p1, &Vx, &Vy, &Vz);
      Tx[p0] += Vx, Ty[p0] += Vy, Tz[p0] += Vz;
      Tx[p1] += Vx, Ty[p1] += Vy, Tz[p1] += Vz;
    }

    
    for(PetscInt p = 0; p < solid->N_total_pts; ++p)  {
      solid->Sp[p] = 0.5*CalcMag( Tx[p], Ty[p], Tz[p]);
      
      solid->Nx[p] =  Ty[p]/(2*solid->Sp[p]);
      solid->Ny[p] = -Tx[p]/(2*solid->Sp[p]);
      solid->Nz[p] = 0.0;
    }
    
    SolidCalculateWeightedCenter(solid, &Gx, &Gy, &Gz);
    for (PetscInt p = 0; p < solid->N_total_pts; ++p) 
    {
      
      CalcVector2(Gx, Gy, Gz, solid->x[p], solid->y[p], solid->z[p], &Vgpx, &Vgpy, &Vgpz);

      
      if ( (Vgpx*solid->Nx[p] + Vgpy*solid->Ny[p]) < 0)
        solid->Nx[p] = -1*solid->Nx[p], solid->Ny[p] = -1*solid->Ny[p], solid->Nz[p] = 0.0;
    }

    PetscFree3(Tx, Ty, Tz); 
  }
  else if (system->Dim == 3 && (solid->Dim == 3 || solid->Dim == 2))  
    
  {
    
    PetscScalar Cx, Cy, Cz; 
    PetscScalar M01x, M01y, M01z, M02x, M02y, M02z, M12x, M12y, M12z; 
    PetscScalar Gx, Gy, Gz; 
    PetscScalar Vgcx, Vgcy, Vgcz; 
    PetscScalar Ntx, Nty, Ntz; 
    PetscInt *N_counter;

    if (solid->N_total_trs <= 0){
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error in SolidCalculateSizeAndNormalAtEachPoint");  }

    PetscMalloc1(solid->N_total_pts, &N_counter); 
    for(PetscInt p = 0; p < solid->N_total_pts; ++p) N_counter[p] = 0;

    SolidCalculateGeometricCenter(solid, &Gx, &Gy, &Gz);

    for(PetscInt t=0; t < solid->N_total_trs; ++t)
    {
      const PetscInt p0 = solid->trs[t][0];  
      const PetscInt p1 = solid->trs[t][1];  
      const PetscInt p2 = solid->trs[t][2];  

      
      CalcElementCentroid(solid->x, solid->y, solid->z, solid->trs[t], 3, &Cx, &Cy, &Cz); 
      
      M01x = (solid->x[p0] + solid->x[p1])/2.0;
      M01y = (solid->y[p0] + solid->y[p1])/2.0;
      M01z = (solid->z[p0] + solid->z[p1])/2.0;
      M02x = (solid->x[p0] + solid->x[p2])/2.0;
      M02y = (solid->y[p0] + solid->y[p2])/2.0;
      M02z = (solid->z[p0] + solid->z[p2])/2.0;
      M12x = (solid->x[p1] + solid->x[p2])/2.0;
      M12y = (solid->y[p1] + solid->y[p2])/2.0;
      M12z = (solid->z[p1] + solid->z[p2])/2.0;

      solid->Sp[p0] += CalcAreaOfTriangle2(Cx, Cy, Cz, M01x, M01y, M01z, solid->x[p0], solid->y[p0], solid->z[p0]);
      solid->Sp[p0] += CalcAreaOfTriangle2(Cx, Cy, Cz, M02x, M02y, M02z, solid->x[p0], solid->y[p0], solid->z[p0]);
      solid->Sp[p1] += CalcAreaOfTriangle2(Cx, Cy, Cz, M01x, M01y, M01z, solid->x[p1], solid->y[p1], solid->z[p1]);
      solid->Sp[p1] += CalcAreaOfTriangle2(Cx, Cy, Cz, M12x, M12y, M12z, solid->x[p1], solid->y[p1], solid->z[p1]);      
      solid->Sp[p2] += CalcAreaOfTriangle2(Cx, Cy, Cz, M02x, M02y, M02z, solid->x[p2], solid->y[p2], solid->z[p2]);
      solid->Sp[p2] += CalcAreaOfTriangle2(Cx, Cy, Cz, M12x, M12y, M12z, solid->x[p2], solid->y[p2], solid->z[p2]);
      
      CalcNormalOfPlane2( solid->x[p0], solid->y[p0], solid->z[p0], solid->x[p1], solid->y[p1], solid->z[p1], 
                          solid->x[p2], solid->y[p2], solid->z[p2], &Ntx, &Nty, &Ntz);
      
      CalcVector2( Gx, Gy, Gz, Cx, Cy, Cz, &Vgcx, &Vgcy,&Vgcz);
      
      if( ((Vgcx*Ntx + Vgcy*Nty + Vgcz*Ntz) < 0))
        Ntx = -1*Ntx, Nty = -1*Nty, Ntz = -1*Ntz;

      solid->Nx[p0] += Ntx; solid->Ny[p0] += Nty;  solid->Nz[p0] += Ntz;
      solid->Nx[p1] += Ntx; solid->Ny[p1] += Nty;  solid->Nz[p1] += Ntz;
      solid->Nx[p2] += Ntx; solid->Ny[p2] += Nty;  solid->Nz[p2] += Ntz;
      N_counter[p0]++, N_counter[p1]++, N_counter[p2]++;
    }

    for(PetscInt p = 0; p < solid->N_total_pts; ++p)
    {
      solid->Nx[p] = solid->Nx[p] / N_counter[p]; 
      solid->Ny[p] = solid->Ny[p] / N_counter[p]; 
      solid->Nz[p] = solid->Nz[p] / N_counter[p];
    }
    PetscFree(N_counter);
  }
  else { 
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error in SolidCalculateSizeAndNormalAtEachPoint");  
  }
  PetscBarrier(NULL); 
  return 0;
}


PetscErrorCode SolidResetVelocity(Solid_Struct *solid) 
{
  for(PetscInt i=0; i<solid->N_total_pts; i++){
    solid->Ux[i] = 0, solid->Uy[i] = 0, solid->Uz[i] = 0;}
  return 0;
}


PetscErrorCode SolidMovePoints(Solid_Struct *solid, const PetscMPIInt rank, const PetscScalar dt,
  PetscInt *ps_mpi, PetscScalar *xs_mpi, PetscScalar *ys_mpi, PetscScalar *zs_mpi )
{

  PetscInt j = 0;
  for(PetscInt p=0; p< solid->N_total_pts; p++)
  {
    if (solid->fluid_elem[p] == -1)
      continue; 
    else
    {
      solid->x[p] += solid->Ux[p]*dt;
      solid->y[p] += solid->Uy[p]*dt;
      solid->z[p] += solid->Uz[p]*dt;
    }
    ps_mpi[j] = p;
    xs_mpi[j] = solid->x[p];
    ys_mpi[j] = solid->y[p];
    zs_mpi[j] = solid->z[p];
    j++;
  }

  if (solid->N_total_pts_loc != j){ 
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Mesh type unknown.");  }

  return 0;
}



PetscErrorCode SolidCalculateVelocityFromDamping(System_Struct *system, Solid_Struct *solid)
{ 
  const PetscScalar c   = solid->forces.c;

  for(PetscInt p=0; p< solid->N_total_pts; p++)      
  {
    
    solid->Ux[p] = solid->forces.Fx[p]/c;
    solid->Uy[p] = solid->forces.Fy[p]/c;
    solid->Uz[p] = solid->forces.Fz[p]/c;
  }
  return 0;
}


PetscErrorCode SolidMovePointsSerial(Solid_Struct *solid, const PetscScalar dt )
{
  
  for(PetscInt p=0; p< solid->N_total_pts; p++)
  {
    solid->x[p] += solid->Ux[p]*dt;
    solid->y[p] += solid->Uy[p]*dt;
    solid->z[p] += solid->Uz[p]*dt;
  }
  PetscPrintf(PETSC_COMM_WORLD,"  Moved all poionts\n");
  return 0;
}


#endif 
