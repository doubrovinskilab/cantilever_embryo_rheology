#ifndef TransferForceToFluid_C
#define TransferForceToFluid_C

PetscErrorCode TransferForceToFluid(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements,
  StokesVariables_Struct *stokesVariables, FullLinearSystem_Stokes_Struct *stokesFullLS,
  Boxes_Struct *boxes, Solid_Struct *solid, BaseFunction *Base_Func)
{
  PetscErrorCode ierr;
  PetscLogDouble Clock1, Clock2;

  PetscTime(&Clock1);
  const PetscInt N    = elements->U_S_inter; 
  const PetscInt Dim  = system->Dim;
  PetscInt  i;
  PetscScalar bases[N];
  PetscScalar Fx_f, Fy_f, Fz_f; 

  PetscInt  N_ins = 0;   
  PetscInt *I_ins;       
  PetscScalar *V_ins;    
  PetscInt j;            

  PetscPrintf(PETSC_COMM_WORLD,"  Transferring forces to fluid on all ranks\n");
  for (PetscInt p = 0; p < solid->N_total_pts; ++p) 
  {
    if (solid->fluid_elem[p] == -1)
      continue; 
    else
      N_ins += Dim*N;
  }

  ierr = PetscMalloc2(N_ins, &I_ins, N_ins, &V_ins); CHKERRQ(ierr);
  for (PetscInt i = 0; i < N_ins; ++i) I_ins[i] = -1, V_ins[i] = 0.0;;

  j = 0;
  for (PetscInt p = 0; p < solid->N_total_pts; ++p) 
  {
    if (solid->fluid_elem[p] == -1)
      continue; 
    else
    {
      const PetscInt    ep   = solid->fluid_elem[p];  
      const PetscScalar Fx_p = solid->forces.Fx[p]; 
      const PetscScalar Fy_p = solid->forces.Fy[p]; 
      const PetscScalar Fz_p = solid->forces.Fz[p]; 

      (*Base_Func)( ep, nodes->x, nodes->y, nodes->z, elements->inter, 
        solid->x[p], solid->y[p], solid->z[p], bases);

      for (i=0; i<N; ++i){  
        Fx_f = Fx_p*bases[i];           
        Fy_f = Fy_p*bases[i];           
        Fz_f = Fz_p*bases[i]*(Dim-2);   

        I_ins[j] = stokesVariables->GI_V[0][elements->inter[ep][i]];    
        V_ins[j] = Fx_f;
        j++;

        I_ins[j] = stokesVariables->GI_V[1][elements->inter[ep][i]];    
        V_ins[j] = Fy_f;
        j++;

        if (Dim>2){
          I_ins[j] = stokesVariables->GI_V[2][elements->inter[ep][i]];  
          V_ins[j] = Fz_f;
          j++;
        }
      }
    }
  }


  if (j!= N_ins){
    SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error in TransferForceToFluid");  } 

  VecSetValues(stokesFullLS->b,  N_ins,    I_ins,  V_ins,    ADD_VALUES);           

  ierr = PetscFree2(I_ins, V_ins); CHKERRQ(ierr);

  ierr = VecAssemblyBegin(stokesFullLS->b); CHKERRQ(ierr);
  ierr = VecAssemblyEnd(stokesFullLS->b); CHKERRQ(ierr);
  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"    Finished transferring forces in %lf seconds.\n", (Clock2 - Clock1) );
  return 0;
}


#endif 
