#ifndef TransferVelocityToSolid_c
#define TransferVelocityToSolid_c


PetscErrorCode SaveBoundaryPointsPosition(Solid_Struct *solid, PetscInt *N_save, PetscInt *p_save, PetscInt *g_save, 
  PetscScalar *x_save, PetscScalar *y_save, PetscScalar *z_save)
{
  PetscInt N_save_tmp =0;
  PetscBool      saved_p= PETSC_FALSE;

  
  for(PetscInt bp=0; bp<solid->N_bound_pts; bp++)
  {
    const PetscInt g = solid->bound_group_pts[bp]; 
    if (  solid->Gr_type[g] == DIST_FIXED || 
          solid->Gr_type[g] == DIST_SPEED || 
          solid->Gr_type[g] == DIST_SPEED_X_LINEAR_Y || 
          solid->Gr_type[g] == AXIS_FIXED || 
          solid->Gr_type[g] == DIST_SPEED_ELLIPSE)
    {
      const PetscInt p = solid->bound_pts[bp];      
      if (solid->fluid_elem[p] == -1)    continue;  

      saved_p = PETSC_FALSE;
      for(PetscInt k=0; k<N_save_tmp; k++) {
        if (p_save[k] == p)  {
          saved_p = PETSC_TRUE; 
          break; }}

      if ( PetscNot(saved_p) ) 
      {
        p_save[N_save_tmp] = p;
        g_save[N_save_tmp] = g;
        x_save[N_save_tmp] = solid->x[p];
        y_save[N_save_tmp] = solid->y[p];
        z_save[N_save_tmp] = solid->z[p];
        N_save_tmp++;
      }
    }
  }

  
  for(PetscInt bl=0; bl<solid->N_bound_lns; bl++)  {
    const PetscInt g = solid->bound_group_lns[bl]; 
    if (  solid->Gr_type[g] == DIST_FIXED || 
          solid->Gr_type[g] == DIST_SPEED || 
          solid->Gr_type[g] == DIST_SPEED_X_LINEAR_Y || 
          solid->Gr_type[g] == AXIS_FIXED || 
          solid->Gr_type[g] == DIST_SPEED_ELLIPSE) 
    {
      const PetscInt l = solid->bound_lns[bl];      
      for(PetscInt i=0; i<2; i++)
      {
        const PetscInt p = solid->lns[l][i];      
        if (solid->fluid_elem[p] == -1)    continue;  

        saved_p = PETSC_FALSE;
        for(PetscInt k=0; k<N_save_tmp; k++) { 
          if (p_save[k] == p)  {            
            saved_p = PETSC_TRUE; 
            break; }}

        if ( PetscNot(saved_p) ) 
        {
          p_save[N_save_tmp] = p;
          g_save[N_save_tmp] = g;
          x_save[N_save_tmp] = solid->x[p];
          y_save[N_save_tmp] = solid->y[p];
          z_save[N_save_tmp] = solid->z[p];
          N_save_tmp++;
        }
      }
    }
  }

  
  for(PetscInt bt=0; bt<solid->N_bound_trs; bt++)  {
    const PetscInt g = solid->bound_group_trs[bt]; 
    if (  solid->Gr_type[g] == DIST_FIXED || 
          solid->Gr_type[g] == DIST_SPEED || 
          solid->Gr_type[g] == DIST_SPEED_X_LINEAR_Y || 
          solid->Gr_type[g] == AXIS_FIXED || 
          solid->Gr_type[g] == DIST_SPEED_ELLIPSE) 
    {
      const PetscInt t = solid->bound_trs[bt];      
      for(PetscInt i=0; i<3; i++)
      {
        const PetscInt p = solid->trs[t][i];      
        if (solid->fluid_elem[p] == -1)    continue;  

        saved_p = PETSC_FALSE;
        for(PetscInt k=0; k<N_save_tmp; k++) { 
          if (p_save[k] == p)  {            
            saved_p = PETSC_TRUE; 
            break; }}

        if ( PetscNot(saved_p) ) 
        {
          p_save[N_save_tmp] = p;
          g_save[N_save_tmp] = g;
          x_save[N_save_tmp] = solid->x[p];
          y_save[N_save_tmp] = solid->y[p];
          z_save[N_save_tmp] = solid->z[p];
          N_save_tmp++;
        }
      }
    }
  }

  *N_save = N_save_tmp;

  return 0;
}

PetscErrorCode AdjustBoundaryPointsSerial(Solid_Struct *solid, const PetscScalar dt, const PetscInt N_save, PetscInt *p_save, PetscInt *g_save, 
  PetscScalar *x_save, PetscScalar *y_save, PetscScalar *z_save)
{ 
  for(PetscInt j=0; j< N_save; j++)      
  {
    const PetscInt     g = g_save[j];
    const PetscInt     p = p_save[j];
    const PetscScalar  x = x_save[j];
    const PetscScalar  y = y_save[j];
    const PetscScalar  z = z_save[j];
    
    if (solid->Gr_type[g]==DIST_FIXED) 
    {
      solid->x[p] = x;
      solid->y[p] = y;
      solid->z[p] = z;
    }
    else if (solid->Gr_type[g]==DIST_SPEED)
    {
      solid->x[p] = x + solid->Gr_value[g][0]*dt; 
      solid->y[p] = y + solid->Gr_value[g][1]*dt; 
      solid->z[p] = z + solid->Gr_value[g][2]*dt; 
    }
    else if (solid->Gr_type[g]==SYMMETRY)
    {
      solid->x[p] = x + solid->Ux[p]*solid->Gr_value[g][0]*dt; 
      solid->y[p] = y + solid->Uy[p]*solid->Gr_value[g][1]*dt; 
      solid->z[p] = z + solid->Uz[p]*solid->Gr_value[g][2]*dt; 
    }
    else if (solid->Gr_type[g]==DIST_SPEED_ELLIPSE)
    {
      const PetscScalar A         = solid->Gr_value[g][0];
      const PetscScalar B         = solid->Gr_value[g][1];
      const PetscScalar DThetaDt  = solid->Gr_value[g][2];
      PetscScalar Theta = CalcEllipsoidAngle(x, y, A, B);
      Theta += dt*DThetaDt;
      solid->x[p] = A*PetscCosScalar(Theta);
      solid->y[p] = B*PetscSinScalar(Theta);
      solid->z[p] = z; 
      
    }
    else if (solid->Gr_type[g]==DIST_SPEED_X_LINEAR_Y)
    {
      solid->x[p] = x + (solid->Gr_value[g][0]*y+solid->Gr_value[g][1])*dt; 
    }
    else if (solid->Gr_type[g]==AXIS_FIXED) 
    {
      const PetscInt a = (PetscInt) solid->Gr_value[g][0];
      if      (a==0)  solid->x[p] = x;
      else if (a==1)  solid->y[p] = y;
      else if (a==2)  solid->z[p] = z;
    }
  }

  
  for (PetscInt i=0; i< solid->forces.BC_FS.N; i++)
  {
    const PetscInt  g = solid->forces.BC_FS.g[i];
    solid->forces.BC_FS.x[i] += solid->Gr_value[g][0]*dt;
    solid->forces.BC_FS.y[i] += solid->Gr_value[g][1]*dt;
    solid->forces.BC_FS.z[i] += solid->Gr_value[g][2]*dt;
  }

  
  for(PetscInt j=0; j< solid->forces.BC_FSE.N; j++)
  {
    const PetscInt g = solid->forces.BC_FSE.g[j];
    const PetscScalar DThetaDt  = solid->Gr_value[g][2];
    solid->forces.BC_FSE.theta[j] += dt*DThetaDt;
  }

  return 0;
}

PetscErrorCode AdjustBoundaryPoints(Solid_Struct *solid, const PetscScalar dt, const PetscInt N_save, PetscInt *p_save, PetscInt *g_save, 
  PetscScalar *x_save, PetscScalar *y_save, PetscScalar *z_save, PetscInt *ps_mpi, PetscScalar *xs_mpi, PetscScalar *ys_mpi, PetscScalar *zs_mpi)
{ 

  AdjustBoundaryPointsSerial( solid, dt, N_save, p_save, g_save, x_save, y_save, z_save);

  PetscBool found_p = PETSC_FALSE;

  for(PetscInt j=0; j< N_save; j++)      
  {
    const PetscInt     p = p_save[j];

    
    PetscInt i = 0;
    found_p = PETSC_FALSE;
    while(PETSC_TRUE) 
    {
      if( ps_mpi[i] == p) {
        found_p = PETSC_TRUE;
        break; }
      i++; 
    }

    if (found_p)
    { xs_mpi[i] = solid->x[p], ys_mpi[i] = solid->y[p], zs_mpi[i] = solid->z[p];}
  }

  return 0;
}


PetscErrorCode BroadCastPoints(System_Struct *system, Solid_Struct *solid,PetscInt *ps_mpi, PetscScalar *xs_mpi, PetscScalar *ys_mpi, PetscScalar *zs_mpi)
{ 
  PetscErrorCode ierr;
  PetscInt        data_size_mpi;            
  PetscInt       *pr_mpi;
  PetscScalar    *xr_mpi, *yr_mpi, *zr_mpi;
  MPI_Request Request;
  MPI_Status  Status;
  PetscPrintf(PETSC_COMM_WORLD,"    Broadcasting solid points location\n");

  PetscBarrier(NULL);

  if (system->size != 1) 
  {

    if (system->rank != 0) 
    { 
      
      data_size_mpi = solid->N_total_pts_loc;

      
      ierr = MPI_Isend(&data_size_mpi, 1, MPIU_INT, 0, system->rank+10000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

      if (data_size_mpi > 0)
      {
        ierr = MPI_Isend(ps_mpi, data_size_mpi, MPIU_INT, 0, system->rank+20000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

        ierr = MPI_Isend(xs_mpi, data_size_mpi, MPIU_SCALAR, 0, system->rank+30000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

        ierr = MPI_Isend(ys_mpi, data_size_mpi, MPIU_SCALAR, 0, system->rank+40000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

        ierr = MPI_Isend(zs_mpi, data_size_mpi, MPIU_SCALAR, 0, system->rank+50000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
      }
    }
    else  
    {
      for(PetscInt r = 1; r<system->size; ++r) 
      {
        data_size_mpi = -1;
        ierr = MPI_Irecv(&data_size_mpi, 1, MPIU_INT,  r, r+10000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

        if (data_size_mpi > 0)
        {
          PetscMalloc4(data_size_mpi, &pr_mpi, data_size_mpi, &xr_mpi, data_size_mpi, &yr_mpi, data_size_mpi, &zr_mpi);
          ierr = MPI_Irecv( pr_mpi, data_size_mpi, MPIU_INT,  r, r+20000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
          ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

          ierr = MPI_Irecv( xr_mpi, data_size_mpi, MPIU_SCALAR,  r, r+30000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
          ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

          ierr = MPI_Irecv( yr_mpi, data_size_mpi, MPIU_SCALAR,  r, r+40000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
          ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

          ierr = MPI_Irecv( zr_mpi, data_size_mpi, MPIU_SCALAR,  r, r+50000, PETSC_COMM_WORLD, &Request); CHKERRMPI(ierr);
          ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

          for (PetscInt j=0;j<data_size_mpi;j++) 
          {
            solid->x[pr_mpi[j]] = xr_mpi[j];
            solid->y[pr_mpi[j]] = yr_mpi[j];
            solid->z[pr_mpi[j]] = zr_mpi[j];
          }
          PetscFree4(pr_mpi, xr_mpi, yr_mpi, zr_mpi);
        }
      }
    }
    PetscBarrier(NULL);

    ierr = MPI_Bcast(solid->x, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    ierr = MPI_Bcast(solid->y, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);
    ierr = MPI_Bcast(solid->z, solid->N_total_pts, MPIU_SCALAR, 0, PETSC_COMM_WORLD ); CHKERRMPI(ierr);

  }

  return 0;
}

PetscErrorCode TransferVelocityFromFluidToSolidPoints(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements,
  StokesVariables_Struct *stokesVariables, FullLinearSystem_Stokes_Struct *stokesFullLS,
  Boxes_Struct *boxes, Solid_Struct *solid, BaseFunction *Base_Func, const PetscScalar dt)
{

  PetscErrorCode ierr;
  PetscLogDouble Clock1, Clock2;
  PetscTime(&Clock1);

  const PetscInt N = elements->U_S_inter;
  PetscScalar *bases;

  
  PetscInt      N_save, N_save_max;
  PetscInt      *p_save, *g_save;
  PetscScalar   *x_save, *y_save, *z_save;

  
  PetscInt    *ps_mpi;
  PetscScalar *xs_mpi, *ys_mpi, *zs_mpi;

  PetscScalar *array_x;
  PetscInt i_low, i_high;
  PetscScalar U[] = {0.0, 0.0, 0.0};

  PetscPrintf(PETSC_COMM_WORLD,"  Transferring velocity from fluid elements to solid-points on all ranks.\n");

  PetscMalloc1(N, &bases);

  ierr = VecGetOwnershipRange(stokesFullLS->x,  &i_low, &i_high); CHKERRQ(ierr); 
  ierr = VecGhostGetLocalForm(stokesFullLS->x,  &stokesFullLS->x_loc); CHKERRQ(ierr);
  ierr = VecGetArray(stokesFullLS->x_loc,       &array_x); CHKERRQ(ierr);

  
  N_save_max = solid->N_bound*3; 
  if (N_save_max>0) 
  {
    PetscMalloc5( N_save_max, &p_save, N_save_max, &g_save, N_save_max, &x_save, N_save_max, &y_save, N_save_max, &z_save);
    SaveBoundaryPointsPosition(solid, &N_save, p_save, g_save, x_save, y_save, z_save);
    if (N_save_max<N_save) {
      PetscErrorPrintf("Error on rank %d. N_save_max %d < N_save.\n", system->rank, N_save_max, N_save);
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Too many boundary points");  } 
  }

  SolidResetVelocity(solid);
  
  for(PetscInt p = 0; p < solid->N_total_pts; ++p) 
  {
    if (solid->fluid_elem[p] == -1)
      continue; 
    else 
    {
      const PetscInt ep = solid->fluid_elem[p]; 

      (*Base_Func)( ep, nodes->x, nodes->y, nodes->z, elements->inter, 
        solid->x[p], solid->y[p], solid->z[p], bases);

      for(PetscInt i=0; i<N; ++i) 
      {      
        const PetscInt n = elements->inter[ep][i]; 
        U[0]=0.0, U[1]=0.0, U[2]=0.0;

        if (nodes->t[n] < 2) 
        {
          for (PetscInt l = 0; l < stokesVariables->Nv_U; l++)  
            U[l] = array_x[stokesVariables->GI_V[l][n] - i_low];

          solid->Ux[p] += bases[i]*U[0];
          solid->Uy[p] += bases[i]*U[1];
          solid->Uz[p] += bases[i]*U[2];
        }
        else 
        {
          PetscInt k = 0;
          for (PetscInt l = 0; l < stokesVariables->Nv_U; l++)
          {
            while(PETSC_TRUE)
            { if(stokesVariables->GI_V[l][n] == stokesVariables->GI_V_ghost[k])  break;
              ++k;}
            U[l]= array_x[k+stokesVariables->N_V_total_loc];
          }

          k = 0;
          while(PETSC_TRUE){
            if (n == nodes->LI_U_ghost[k]) break;
              ++k;}

          solid->Ux[p] += bases[i]*U[0];
          solid->Uy[p] += bases[i]*U[1];
          solid->Uz[p] += bases[i]*U[2];
        }
      }
    }
  }
  
  ierr = VecRestoreArray(stokesFullLS->x_loc, &array_x); CHKERRQ(ierr);
  ierr = VecGhostRestoreLocalForm(stokesFullLS->x, &stokesFullLS->x_loc); CHKERRQ(ierr);

  PetscFree( bases);

  
  PetscMalloc4(solid->N_total_pts_loc, &ps_mpi, solid->N_total_pts_loc, &xs_mpi, solid->N_total_pts_loc, &ys_mpi, solid->N_total_pts_loc, &zs_mpi);
  SolidMovePoints(solid, system->rank, dt, ps_mpi, xs_mpi, ys_mpi, zs_mpi);

  
  if(N_save_max > 0) {
    AdjustBoundaryPoints(solid, dt, N_save, p_save, g_save, x_save, y_save, z_save, ps_mpi, xs_mpi, ys_mpi, zs_mpi);
    PetscFree5(p_save, g_save, x_save, y_save, z_save);}

  
  ConstrainCorrectPoints(system, solid, ps_mpi, xs_mpi, ys_mpi, zs_mpi);

  BroadCastPoints( system, solid, ps_mpi, xs_mpi, ys_mpi, zs_mpi); 
  PetscFree4(ps_mpi, xs_mpi, ys_mpi, zs_mpi);

  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"    Finished transferring velocities in %lf seconds.\n", (Clock2 - Clock1) );
  return 0;
}

#endif 
