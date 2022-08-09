#ifndef ElementalMatrixBoundaryPressureNeumann_C
#define ElementalMatrixBoundaryPressureNeumann_C


PetscErrorCode Calc2DP1_Pnlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, 
  PetscInt **e_inter,const PetscInt N_inter_loc, PetscScalar **Pnlp, const PetscScalar *DC, const PetscMPIInt rank)
{
  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  PetscInt n2;      
  PetscBool Found = PETSC_FALSE;

  for (PetscInt e=0; e<N_inter_loc; e++)
  {
    const PetscInt ni0 = e_inter[e][0];
    const PetscInt ni1 = e_inter[e][1];
    const PetscInt ni2 = e_inter[e][2];
    if ((ni0 == n0 && ni1 == n1) || (ni0 == n1 && ni1 == n0)) {
      n2 = ni2; Found = PETSC_TRUE; break; }
    if ((ni0 == n0 && ni2 == n1) || (ni0 == n1 && ni2 == n0)) {
      n2 = ni1; Found = PETSC_TRUE; break; }
    if ((ni1 == n0 && ni2 == n1) || (ni1 == n1 && ni2 == n0)) {
      n2 = ni0; Found = PETSC_TRUE; break; }
  }
  
  if (PetscNot(Found))
  { PetscErrorPrintf("Error in Calc2DP1_Pnlp. 'n2' could not be found for boundary nodes n0=%d, n1=%d, n2=%d at rank %d\n", n0, n1, rank);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Error in Calc2DP1_Pnlp"); }

  
  const PetscScalar x0 = n_x[n0], y0 = n_y[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2];
  
  const PetscScalar Delta = x0*y1-x0*y2-x1*y0+x1*y2+x2*y0-x2*y1;
  
  PetscScalar Nx, Ny, Nmag; 
  PetscScalar nx, ny;       

  const PetscScalar Dx = x1-x0,  Dy = y1-y0;   
  const PetscScalar L = PetscSqrtScalar(Dx*Dx + Dy*Dy); 

  if (Dx != 0.0)  Nx = -Dy/Dx, Ny = 1.0; 
  else            Nx = -1.0, Ny = Dx/Dy; 
  Nmag = PetscSqrtScalar(Nx*Nx + Ny*Ny);
  
  nx = Nx / Nmag;
  ny = Ny / Nmag;
  
  const PetscScalar Cx = x2 - DC[0], Cy = y2 - DC[1]; 
  const PetscScalar DOT = Cx*nx + Cy*ny;              

  if (DOT < 0) nx = -1*nx, ny = -1*ny;
 
  
  Pnlp[0][0] =    L*(y1 - y2)/(2*Delta)*nx;
  Pnlp[0][1] =   -L*(y0 - y2)/(2*Delta)*nx;
  Pnlp[1][0] =    L*(y1 - y2)/(2*Delta)*nx;
  Pnlp[1][1] =   -L*(y0 - y2)/(2*Delta)*nx;

  
  Pnlp[0][0] +=  -L*(x1 - x2)/(2*Delta)*ny;
  Pnlp[0][1] +=   L*(x0 - x2)/(2*Delta)*ny;
  Pnlp[1][0] +=  -L*(x1 - x2)/(2*Delta)*ny;
  Pnlp[1][1] +=   L*(x0 - x2)/(2*Delta)*ny;

  return 0;
}


PetscErrorCode Calc3DP1_Pnlp(const PetscInt en, PetscScalar *n_x, PetscScalar *n_y, PetscScalar *n_z,  PetscInt **e_bound, 
  PetscInt **e_inter,const PetscInt N_inter_loc, PetscScalar **Pnlp, const PetscScalar *DC, const PetscMPIInt rank)
{

  const PetscInt n0 = e_bound[en][0];
  const PetscInt n1 = e_bound[en][1];
  const PetscInt n2 = e_bound[en][2];
  PetscInt n3;
  PetscBool Found = PETSC_FALSE;

  
  for (PetscInt e=0; e<N_inter_loc; e++)
  {
    const PetscInt ni0 = e_inter[e][0];
    const PetscInt ni1 = e_inter[e][1];
    const PetscInt ni2 = e_inter[e][2];
    const PetscInt ni3 = e_inter[e][3];
    
    

    if ( (ni0 == n0 && ni1 == n1 && ni2 == n2) || 
         (ni0 == n1 && ni1 == n0 && ni2 == n2) || 
         (ni0 == n2 && ni1 == n0 && ni2 == n1) || 
         (ni0 == n2 && ni1 == n1 && ni2 == n0) || 
         (ni0 == n0 && ni1 == n2 && ni2 == n1) || 
         (ni0 == n1 && ni1 == n2 && ni2 == n0) ) {
          n3 = ni3; Found = PETSC_TRUE; break; }
    if ( (ni0 == n0 && ni1 == n1 && ni3 == n2) || 
         (ni0 == n1 && ni1 == n0 && ni3 == n2) || 
         (ni0 == n2 && ni1 == n0 && ni3 == n1) || 
         (ni0 == n2 && ni1 == n1 && ni3 == n0) || 
         (ni0 == n0 && ni1 == n2 && ni3 == n1) || 
         (ni0 == n1 && ni1 == n2 && ni3 == n0) ) {
          n3 = ni2; Found = PETSC_TRUE; break; }
    if ( (ni0 == n0 && ni3 == n1 && ni2 == n2) || 
         (ni0 == n1 && ni3 == n0 && ni2 == n2) || 
         (ni0 == n2 && ni3 == n0 && ni2 == n1) || 
         (ni0 == n2 && ni3 == n1 && ni2 == n0) || 
         (ni0 == n0 && ni3 == n2 && ni2 == n1) || 
         (ni0 == n1 && ni3 == n2 && ni2 == n0) ) {
          n3 = ni1; Found = PETSC_TRUE; break; }
    if ( (ni3 == n0 && ni1 == n1 && ni2 == n2) || 
         (ni3 == n1 && ni1 == n0 && ni2 == n2) || 
         (ni3 == n2 && ni1 == n0 && ni2 == n1) || 
         (ni3 == n2 && ni1 == n1 && ni2 == n0) || 
         (ni3 == n0 && ni1 == n2 && ni2 == n1) || 
         (ni3 == n1 && ni1 == n2 && ni2 == n0) ) {
          n3 = ni0; Found = PETSC_TRUE; break; }
  }

  if (PetscNot(Found))
  { PetscErrorPrintf("Error in Calc3DP1_Pnlp. 'n3' could not be found for boundary nodes n0=%d, n1=%d, n2=%d at rank %d\n", n0, n1, n2, rank);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_SUP, "Error in Calc3DP1_Pnlp"); }
  
  const PetscScalar x0 = n_x[n0], y0 = n_y[n0], z0 = n_z[n0];
  const PetscScalar x1 = n_x[n1], y1 = n_y[n1], z1 = n_z[n1];
  const PetscScalar x2 = n_x[n2], y2 = n_y[n2], z2 = n_z[n2];
  const PetscScalar x3 = n_x[n3], y3 = n_y[n3], z3 = n_z[n3];

  const PetscScalar Delta = -x0*y1*z2 + x0*y1*z3 + x0*y2*z1 - x0*y2*z3 - x0*y3*z1 + x0*y3*z2 + x1*y0*z2 
     -x1*y0*z3 - x1*y2*z0 + x1*y2*z3 + x1*y3*z0 - x1*y3*z2 - x2*y0*z1 + x2*y0*z3  
     +x2*y1*z0 - x2*y1*z3 - x2*y3*z0 + x2*y3*z1 + x3*y0*z1 - x3*y0*z2 - x3*y1*z0  
     +x3*y1*z2 + x3*y2*z0 - x3*y2*z1;

  
  PetscScalar Nx, Ny, Nz, Nmag; 
  PetscScalar nx, ny, nz;       

  const PetscScalar xc = (x0 + x1 + x2)/3.0;
  const PetscScalar yc = (y0 + y1 + y2)/3.0;
  const PetscScalar zc = (z0 + z1 + z2)/3.0;
  
  const PetscScalar Dx10 = x1-x0, Dy10 = y1-y0, Dz10 = z1-z0;
  const PetscScalar Dx20 = x2-x0, Dy20 = y2-y0, Dz20 = z2-z0;

  const PetscScalar Area = 0.5*PetscSqrtScalar( PetscPowScalar((Dx10*Dy20 - Dx20*Dy10),2) 
    + PetscPowScalar((-Dx10*Dz20 + Dx20*Dz10),2) + PetscPowScalar((Dy10*Dz20 - Dy20*Dz10),2) );

  const PetscScalar C0yz =  ( y1*z2 - y1*z3 - y2*z1 + y2*z3 + y3*z1 - y3*z2);
  const PetscScalar C1yz =  (-y0*z2 + y0*z3 + y2*z0 - y2*z3 - y3*z0 + y3*z2);
  const PetscScalar C2yz = -(-y0*z1 + y0*z3 + y1*z0 - y1*z3 - y3*z0 + y3*z1);

  const PetscScalar C0xz = -(x1*z2 - x1*z3 - x2*z1 + x2*z3 + x3*z1 - x3*z2);
  const PetscScalar C1xz =  (x0*z2 - x0*z3 - x2*z0 + x2*z3 + x3*z0 - x3*z2);
  const PetscScalar C2xz = -(x0*z1 - x0*z3 - x1*z0 + x1*z3 + x3*z0 - x3*z1);

  const PetscScalar C0xy =  (x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2);
  const PetscScalar C1xy = -(x0*y2 - x0*y3 - x2*y0 + x2*y3 + x3*y0 - x3*y2);
  const PetscScalar C2xy =  (x0*y1 - x0*y3 - x1*y0 + x1*y3 + x3*y0 - x3*y1);

  
  Nx =  Dy10*Dz20 - Dy20*Dz10;
  Ny = -Dx10*Dz20 + Dx20*Dz10;
  Nz =  Dx10*Dy20 - Dx20*Dy10;
  Nmag = PetscSqrtScalar(Nx*Nx + Ny*Ny + Nz*Nz); 
  nx = Nx / Nmag;
  ny = Ny / Nmag;
  nz = Nz / Nmag;

  
  const PetscScalar Cx = xc - DC[0], Cy = yc - DC[1], Cz = zc - DC[2]; 
  const PetscScalar DOT = Cx*nx + Cy*ny + Cz*nz;                       

  if (DOT < 0) nx = -1*nx, ny = -1*ny, nz=-1*nz;

  
  
  Pnlp[0][0] =   -Area*C0yz/(3*Delta)*nx;
  Pnlp[0][1] =   -Area*C0yz/(3*Delta)*nx;
  Pnlp[0][2] =   -Area*C0yz/(3*Delta)*nx; 
  Pnlp[1][0] =   -Area*C1yz/(3*Delta)*nx; 
  Pnlp[1][1] =   -Area*C1yz/(3*Delta)*nx; 
  Pnlp[1][2] =   -Area*C1yz/(3*Delta)*nx;  
  Pnlp[2][0] =   -Area*C2yz/(3*Delta)*nx; 
  Pnlp[2][1] =   -Area*C2yz/(3*Delta)*nx;      
  Pnlp[2][2] =   -Area*C2yz/(3*Delta)*nx;  

  
  Pnlp[0][0] +=  -Area*C0xz/(3*Delta)*ny;
  Pnlp[0][1] +=  -Area*C0xz/(3*Delta)*ny;
  Pnlp[0][2] +=  -Area*C0xz/(3*Delta)*ny;
  Pnlp[1][0] +=  -Area*C1xz/(3*Delta)*ny; 
  Pnlp[1][1] +=  -Area*C1xz/(3*Delta)*ny;
  Pnlp[1][2] +=  -Area*C1xz/(3*Delta)*ny;
  Pnlp[2][0] +=  -Area*C2xz/(3*Delta)*ny;
  Pnlp[2][1] +=  -Area*C2xz/(3*Delta)*ny; 
  Pnlp[2][2] +=  -Area*C2xz/(3*Delta)*ny;   

  
  Pnlp[0][0] +=  -Area*C0xy/(3*Delta)*nz;  
  Pnlp[0][1] +=  -Area*C0xy/(3*Delta)*nz; 
  Pnlp[0][2] +=  -Area*C0xy/(3*Delta)*nz;    
  Pnlp[1][0] +=  -Area*C1xy/(3*Delta)*nz;  
  Pnlp[1][1] +=  -Area*C1xy/(3*Delta)*nz;
  Pnlp[1][2] +=  -Area*C1xy/(3*Delta)*nz;
  Pnlp[2][0] +=  -Area*C2xy/(3*Delta)*nz; 
  Pnlp[2][1] +=  -Area*C2xy/(3*Delta)*nz; 
  Pnlp[2][2] +=  -Area*C2xy/(3*Delta)*nz;  

  return 0;
}

#endif