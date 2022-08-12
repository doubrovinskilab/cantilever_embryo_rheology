#ifndef ElementalMatrixBoundaryPressure_C
#define ElementalMatrixBoundaryPressure_C

PetscErrorCode Calc2DP2P1_Plp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, const PetscInt v, 
  PetscScalar **Plp, const PetscScalar *DC)
{

  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscInt n2 = e_bound[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];
  
  const PetscScalar Dx = x1-x0,  Dy = y1-y0;   
  PetscScalar Nx, Ny;

  if (Dx != 0.0)  Nx = -Dy/Dx, Ny = 1.0; 
  else            Nx = -1.0, Ny = Dx/Dy; 

  if (v==0) 
  {
    Plp[0][0] = -Dy/6.;  Plp[0][1] =     0.;
    Plp[1][0] =     0.;  Plp[1][1] = -Dy/6.;
    Plp[2][0] = -Dy/3.;  Plp[2][1] = -Dy/3.;
  }
  else if (v==1) 
  {
    Plp[0][0] =  Dx/6.;  Plp[0][1] =     0.;
    Plp[1][0] =     0.;  Plp[1][1] =  Dx/6.;
    Plp[2][0] =  Dx/3.;  Plp[2][1] =  Dx/3.;
  }

  const PetscScalar Cx = x2 - DC[0], Cy = y2 - DC[1]; 
  const PetscScalar DOT = Cx*Nx + Cy*Ny;              

  if (DOT < 0) 
  {
    Plp[0][0] = -1*Plp[0][0]; Plp[0][1] = -1*Plp[0][1];
    Plp[1][0] = -1*Plp[1][0]; Plp[1][1] = -1*Plp[1][1];
    Plp[2][0] = -1*Plp[2][0]; Plp[2][1] = -1*Plp[2][1];
  }

  return 0;
}

PetscErrorCode Calc2DP1bP1_Plp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, const PetscInt v, 
  PetscScalar **Plp, const PetscScalar *DC)
{

  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  
  const PetscScalar Dx = x1-x0,  Dy = y1-y0;   
  PetscScalar Nx, Ny;

  if (Dx!=0.0)  Nx = -Dy/Dx,  Ny = 1.0;   
  else          Nx = -1.0,    Ny = Dx/Dy; 

  if (v==0) 
  {
    Plp[0][0] = -Dy/3.;  Plp[0][1] = -Dy/6.;
    Plp[1][0] = -Dy/6.;  Plp[1][1] = -Dy/3.;
  }
  else if (v==1) 
  {
    Plp[0][0] =  Dx/3.;  Plp[0][1] =  Dx/6.; 
    Plp[1][0] =  Dx/6.;  Plp[1][1] =  Dx/3.;
  }

  
  const PetscScalar Cx = (x0+x1)/2.0 - DC[0], Cy = (y0+y1)/2.0 - DC[1]; 
  const PetscScalar DOT = Cx*Nx + Cy*Ny;              

  if (DOT < 0) 
  {
    Plp[0][0] = -1.*Plp[0][0],     Plp[0][1] = -1.*Plp[0][1];
    Plp[1][0] = -1.*Plp[1][0],     Plp[1][1] = -1.*Plp[1][1];
  }
  return 0;
}


PetscErrorCode Calc3DP2P1_Plp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, const PetscInt v, 
  PetscScalar **Plp, const PetscScalar *DC)
{

  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscInt n2 = e_bound[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar xc = (x0 + x1 + x2)/3.0;
  const PetscScalar yc = (y0 + y1 + y2)/3.0;
  const PetscScalar zc = (z0 + z1 + z2)/3.0;
  
  const PetscScalar Dx10 = x1-x0, Dy10 = y1-y0, Dz10 = z1-z0;
  const PetscScalar Dx20 = x2-x0, Dy20 = y2-y0, Dz20 = z2-z0;

  const PetscScalar Area = 0.5*PetscSqrtScalar( PetscPowScalar((Dx10*Dy20 - Dx20*Dy10),2) + PetscPowScalar((-Dx10*Dz20 + Dx20*Dz10),2) + PetscPowScalar((Dy10*Dz20 - Dy20*Dz10),2) );

  const PetscScalar C3yz = -y0*z1 + y0*z2 + y1*z0 - y1*z2 - y2*z0 + y2*z1;
  const PetscScalar C3xz =  x0*z1 - x0*z2 - x1*z0 + x1*z2 + x2*z0 - x2*z1;
  const PetscScalar C3xy = -x0*y1 + x0*y2 + x1*y0 - x1*y2 - x2*y0 + x2*y1;

  
  PetscScalar Nx =  Dy10*Dz20 - Dy20*Dz10;
  PetscScalar Ny = -Dx10*Dz20 + Dx20*Dz10;
  PetscScalar Nz =  Dx10*Dy20 - Dx20*Dy10;

  
  const PetscScalar Cx = xc - DC[0], Cy = yc - DC[1], Cz = zc - DC[2]; 
  const PetscScalar DOT = Cx*Nx + Cy*Ny + Cz*Nz;              

  PetscScalar K = 1.0; 
  if (DOT < 0) K = -1.0;

  if (v==0) 
  {
    Plp[0][0] =   -K*Area*C3yz/30.;
    Plp[0][1] =    K*Area*C3yz/60.;
    Plp[0][2] =    K*Area*C3yz/60.;
    Plp[1][0] =    K*Area*C3yz/60.;
    Plp[1][1] =   -K*Area*C3yz/30.;
    Plp[1][2] =    K*Area*C3yz/60.;
    Plp[2][0] =    K*Area*C3yz/60.;
    Plp[2][1] =    K*Area*C3yz/60.;
    Plp[2][2] =   -K*Area*C3yz/30.;
    Plp[3][0] = -2*K*Area*C3yz/15.;
    Plp[3][1] = -2*K*Area*C3yz/15.;
    Plp[3][2] =   -K*Area*C3yz/15.;
    Plp[4][0] =   -K*Area*C3yz/15.;
    Plp[4][1] = -2*K*Area*C3yz/15.;
    Plp[4][2] = -2*K*Area*C3yz/15.;
    Plp[5][0] = -2*K*Area*C3yz/15.;
    Plp[5][1] =   -K*Area*C3yz/15.;
    Plp[5][2] = -2*K*Area*C3yz/15.;
  }
  else if (v==1) 
  {
    Plp[0][0] =   -K*Area*C3xz/30.;
    Plp[0][1] =    K*Area*C3xz/60.;
    Plp[0][2] =    K*Area*C3xz/60.;
    Plp[1][0] =    K*Area*C3xz/60.;
    Plp[1][1] =   -K*Area*C3xz/30.;
    Plp[1][2] =    K*Area*C3xz/60.;
    Plp[2][0] =    K*Area*C3xz/60.;
    Plp[2][1] =    K*Area*C3xz/60.;
    Plp[2][2] =   -K*Area*C3xz/30.;
    Plp[3][0] = -2*K*Area*C3xz/15.;
    Plp[3][1] = -2*K*Area*C3xz/15.;
    Plp[3][2] =   -K*Area*C3xz/15.;
    Plp[4][0] =   -K*Area*C3xz/15.;
    Plp[4][1] = -2*K*Area*C3xz/15.;
    Plp[4][2] = -2*K*Area*C3xz/15.;
    Plp[5][0] = -2*K*Area*C3xz/15.;
    Plp[5][1] =   -K*Area*C3xz/15.;
    Plp[5][2] = -2*K*Area*C3xz/15.;
  }
  else if (v==2) 
  {
    Plp[0][0] =   -K*Area*C3xy/30.;
    Plp[0][1] =    K*Area*C3xy/60.;
    Plp[0][2] =    K*Area*C3xy/60.;
    Plp[1][0] =    K*Area*C3xy/60.;
    Plp[1][1] =   -K*Area*C3xy/30.;
    Plp[1][2] =    K*Area*C3xy/60.;
    Plp[2][0] =    K*Area*C3xy/60.;
    Plp[2][1] =    K*Area*C3xy/60.;
    Plp[2][2] =   -K*Area*C3xy/30.;
    Plp[3][0] = -2*K*Area*C3xy/15.;
    Plp[3][1] = -2*K*Area*C3xy/15.;
    Plp[3][2] =   -K*Area*C3xy/15.;
    Plp[4][0] =   -K*Area*C3xy/15.;
    Plp[4][1] = -2*K*Area*C3xy/15.;
    Plp[4][2] = -2*K*Area*C3xy/15.;
    Plp[5][0] = -2*K*Area*C3xy/15.;
    Plp[5][1] =   -K*Area*C3xy/15.;
    Plp[5][2] = -2*K*Area*C3xy/15.;
  }
  return 0;
}

PetscErrorCode Calc3DP1bP1_Plp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, const PetscInt v, 
  PetscScalar **Plp, const PetscScalar *DC)
{

  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscInt n2 = e_bound[en][2];

  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar xc = (x0 + x1 + x2)/3.0;
  const PetscScalar yc = (y0 + y1 + y2)/3.0;
  const PetscScalar zc = (z0 + z1 + z2)/3.0;
  
  const PetscScalar Dx10 = x1-x0, Dy10 = y1-y0, Dz10 = z1-z0;
  const PetscScalar Dx20 = x2-x0, Dy20 = y2-y0, Dz20 = z2-z0;

  const PetscScalar Area = 0.5*PetscSqrtScalar( PetscPowScalar((Dx10*Dy20 - Dx20*Dy10),2) + PetscPowScalar((-Dx10*Dz20 + Dx20*Dz10),2) + PetscPowScalar((Dy10*Dz20 - Dy20*Dz10),2) );

  const PetscScalar C3yz = -y0*z1 + y0*z2 + y1*z0 - y1*z2 - y2*z0 + y2*z1;
  const PetscScalar C3xz =  x0*z1 - x0*z2 - x1*z0 + x1*z2 + x2*z0 - x2*z1;
  const PetscScalar C3xy = -x0*y1 + x0*y2 + x1*y0 - x1*y2 - x2*y0 + x2*y1;

  
  PetscScalar Nx =  Dy10*Dz20 - Dy20*Dz10;
  PetscScalar Ny = -Dx10*Dz20 + Dx20*Dz10;
  PetscScalar Nz =  Dx10*Dy20 - Dx20*Dy10;

  
  const PetscScalar Cx = xc - DC[0], Cy = yc - DC[1], Cz = zc - DC[2];  
  const PetscScalar DOT = Cx*Nx + Cy*Ny + Cz*Nz;                        

  PetscScalar K = 1.0; 
  if (DOT < 0) K = -1.0;

  if (v==0) 
  {
    Plp[0][0] = -K*Area*C3yz/6.;
    Plp[0][1] = -K*Area*C3yz/12.;
    Plp[0][2] = -K*Area*C3yz/12.;
    Plp[1][0] = -K*Area*C3yz/12.;
    Plp[1][1] = -K*Area*C3yz/6.;
    Plp[1][2] = -K*Area*C3yz/12.;
    Plp[2][0] = -K*Area*C3yz/12.;
    Plp[2][1] = -K*Area*C3yz/12.;
    Plp[2][2] = -K*Area*C3yz/6.;
  }
  else if (v==1) 
  {
    Plp[0][0] = -K*Area*C3xz/6.;
    Plp[0][1] = -K*Area*C3xz/12.;
    Plp[0][2] = -K*Area*C3xz/12.;
    Plp[1][0] = -K*Area*C3xz/12.;
    Plp[1][1] = -K*Area*C3xz/6.;
    Plp[1][2] = -K*Area*C3xz/12.;
    Plp[2][0] = -K*Area*C3xz/12.;
    Plp[2][1] = -K*Area*C3xz/12.;
    Plp[2][2] = -K*Area*C3xz/6.;
  }
  else if (v==2) 
  {
    Plp[0][0] = -K*Area*C3xy/6.;
    Plp[0][1] = -K*Area*C3xy/12.;
    Plp[0][2] = -K*Area*C3xy/12.;
    Plp[1][0] = -K*Area*C3xy/12.;
    Plp[1][1] = -K*Area*C3xy/6.;
    Plp[1][2] = -K*Area*C3xy/12.;
    Plp[2][0] = -K*Area*C3xy/12.;
    Plp[2][1] = -K*Area*C3xy/12.;
    Plp[2][2] = -K*Area*C3xy/6.;
  }
  return 0;
}
#endif