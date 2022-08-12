#ifndef Math_Functions_C
#define Math_Functions_C

PetscErrorCode CalcElementCentroid(PetscScalar *x, PetscScalar *y, PetscScalar *z, PetscInt *elem, PetscInt elem_size, 
  PetscScalar *Ex, PetscScalar *Ey, PetscScalar *Ez) 
{
  *Ex = 0, *Ey = 0, *Ez = 0;
  for (PetscInt j = 0; j < elem_size; j++) 
  {
    const PetscInt n = elem[j];
    *Ex += x[n];
    *Ey += y[n];
    *Ez += z[n];
  }
  *Ex = *Ex/elem_size;
  *Ey = *Ey/elem_size;
  *Ez = *Ez/elem_size;

  return 0;
}

PetscScalar CalcDistance(PetscScalar *x, PetscScalar *y, PetscScalar *z, const PetscInt p0, const PetscInt p1) 
{
 PetscScalar Dx = x[p1]-x[p0];
 PetscScalar Dy = y[p1]-y[p0];
 PetscScalar Dz = z[p1]-z[p0];
 return PetscSqrtScalar( Dx*Dx + Dy*Dy + Dz*Dz );
}


PetscScalar CalcDistance2(const PetscScalar p0x, const PetscScalar p0y, const PetscScalar p0z, 
  const PetscScalar p1x, const PetscScalar p1y, const PetscScalar p1z) 
{
 PetscScalar Dx = p1x-p0x;
 PetscScalar Dy = p1y-p0y;
 PetscScalar Dz = p1z-p0z;
 return PetscSqrtScalar( Dx*Dx + Dy*Dy + Dz*Dz );
}

PetscScalar CalcEllipsoidAngle(const PetscScalar x, const PetscScalar y, const PetscScalar A, const PetscScalar B)
{
  const PetscScalar   PI = 3.14159265358979323846; 
  PetscScalar Theta;
  if (x==0) {
    if      (y<0)   Theta = -PI/2;
    else if (y>0)   Theta =  PI/2;
    else if (y==0)  Theta =   0.0; }
  else {
    Theta = PetscAtanScalar((y*A)/(x*B));
    if (x<0) {
      if      (y<0)  Theta = Theta - PI;  
      else if (y>0)  Theta = Theta + PI;
      else if (y==0) Theta = PI;  }
  }
  return Theta;
}



void CalcVector( PetscScalar *x, PetscScalar *y, PetscScalar *z, const PetscInt p0, const PetscInt p1,
  PetscScalar *Vx, PetscScalar *Vy, PetscScalar *Vz)
{
  *Vx = x[p1] - x[p0];
  *Vy = y[p1] - y[p0];
  *Vz = z[p1] - z[p0];
}

void CalcVector2(const PetscScalar p0x, const PetscScalar p0y, const PetscScalar p0z,
  const PetscScalar p1x, const PetscScalar p1y, const PetscScalar p1z, PetscScalar *Vx, PetscScalar *Vy, PetscScalar *Vz)
{
  *Vx = p1x - p0x;
  *Vy = p1y - p0y;
  *Vz = p1z - p0z;
}

PetscScalar CalcMag( const PetscScalar Vx, const PetscScalar Vy, const PetscScalar Vz){
  return PetscSqrtScalar( Vx*Vx + Vy*Vy + Vz*Vz); }

PetscScalar CalcDotProduct( const PetscScalar V1x, const PetscScalar V1y, const PetscScalar V1z,
                            const PetscScalar V2x, const PetscScalar V2y, const PetscScalar V2z){
  return V1x*V2x + V1y*V2y + V1z*V2z;
}


void CalcCrossProduct( const PetscScalar V1x, const PetscScalar V1y, const PetscScalar V1z,
                       const PetscScalar V2x, const PetscScalar V2y, const PetscScalar V2z,
                             PetscScalar *Vx,       PetscScalar *Vy,       PetscScalar *Vz)
{
  *Vx =   V1y*V2z - V2y*V1z;
  *Vy = -(V1x*V2z - V2x*V1z);
  *Vz =   V1x*V2y - V2x*V1y;
}

void CalcNormalVector( PetscScalar *x, PetscScalar *y, PetscScalar *z, const PetscInt p0, const PetscInt p1,
  PetscScalar *Vx, PetscScalar *Vy, PetscScalar *Vz)
{
  *Vx = x[p1] - x[p0];
  *Vy = y[p1] - y[p0];
  *Vz = z[p1] - z[p0];
  PetscScalar Vm = CalcMag(*Vx, *Vy, *Vz);
  *Vx = *Vx / Vm;
  *Vy = *Vy / Vm;
  *Vz = *Vz / Vm;
}

void CalcNormalVector2( const PetscScalar p0x, const PetscScalar p0y, const PetscScalar p0z,
  const PetscScalar p1x, const PetscScalar p1y, const PetscScalar p1z, PetscScalar *Vx, PetscScalar *Vy, PetscScalar *Vz)
{
  *Vx = p1x - p0x;
  *Vy = p1y - p0y;
  *Vz = p1z - p0z;
  PetscScalar Vm = CalcMag(*Vx, *Vy, *Vz);
  *Vx = *Vx / Vm;
  *Vy = *Vy / Vm;
  *Vz = *Vz / Vm;
}

void CalcNormalOfPlane( PetscScalar *x,   PetscScalar *y, PetscScalar *z, const PetscInt p0, const PetscInt p1, const PetscInt p2,
                        PetscScalar *Ntx, PetscScalar *Nty, PetscScalar *Ntz)
{
  PetscScalar V1x, V1y, V1z;
  PetscScalar V2x, V2y, V2z;
  CalcVector( x, y, z, p0, p1, &V1x, &V1y, &V1z);
  CalcVector( x, y, z, p0, p2, &V2x, &V2y, &V2z);
  CalcCrossProduct( V1x, V1y, V1z, V2x, V2y, V2z, Ntx, Nty, Ntz);
  PetscScalar Ntm = CalcMag( *Ntx, *Nty, *Ntz);
  *Ntx = *Ntx / Ntm;
  *Nty = *Nty / Ntm;
  *Ntz = *Ntz / Ntm;
}

void CalcNormalOfPlane2(const PetscInt p0x, const PetscInt p0y, const PetscInt p0z,
                        const PetscInt p1x, const PetscInt p1y, const PetscInt p1z,
                        const PetscInt p2x, const PetscInt p2y, const PetscInt p2z,
                        PetscScalar *Ntx, PetscScalar *Nty, PetscScalar *Ntz)
{
  PetscScalar V1x, V1y, V1z;
  PetscScalar V2x, V2y, V2z;
  CalcVector2(p0x, p0y, p0z, p1x, p1y, p1z, &V1x, &V1y, &V1z);
  CalcVector2(p0x, p0y, p0z, p2x, p2y, p2z, &V2x, &V2y, &V2z);
  CalcCrossProduct( V1x, V1y, V1z, V2x, V2y, V2z, Ntx, Nty, Ntz);
  PetscScalar Ntm = CalcMag( *Ntx, *Nty, *Ntz);
  *Ntx = *Ntx / Ntm;
  *Nty = *Nty / Ntm;
  *Ntz = *Ntz / Ntm;
}

PetscScalar CalcAreaOfTriangle(PetscScalar *x, PetscScalar *y, PetscScalar *z, const PetscInt p0, const PetscInt p1, const PetscInt p2)
{
  PetscScalar V1x, V1y, V1z;
  PetscScalar V2x, V2y, V2z;
  PetscScalar  Vx,  Vy,  Vz;
  CalcVector( x, y, z, p0, p1, &V1x, &V1y, &V1z);
  CalcVector( x, y, z, p0, p2, &V2x, &V2y, &V2z);
  CalcCrossProduct( V1x, V1y, V1z, V2x, V2y, V2z, &Vx, &Vy, &Vz);
  return 0.5*CalcMag(Vx, Vy, Vz);
}

PetscScalar CalcAreaOfTriangle2(const PetscScalar p0x, const PetscScalar p0y, const PetscScalar p0z, 
  const PetscScalar p1x, const PetscScalar p1y, const PetscScalar p1z, 
  const PetscScalar p2x, const PetscScalar p2y, const PetscScalar p2z)
{
  PetscScalar V1x, V1y, V1z;
  PetscScalar V2x, V2y, V2z;
  PetscScalar  Vx,  Vy,  Vz;
  CalcVector2(p0x, p0y, p0z, p1x, p1y, p1z, &V1x, &V1y, &V1z);
  CalcVector2(p0x, p0y, p0z, p2x, p2y, p2z, &V2x, &V2y, &V2z);
  CalcCrossProduct( V1x, V1y, V1z, V2x, V2y, V2z, &Vx, &Vy, &Vz);
  return 0.5*CalcMag(Vx, Vy, Vz);
}

PetscScalar CalcEllipseEquation(const PetscScalar x, const PetscScalar y, const PetscScalar z, const PetscScalar A, const PetscScalar B, const PetscScalar C){
  return (x*x)/(A*A) + (y*y)/(B*B) + (z*z)/(C*C);
}


#endif 
