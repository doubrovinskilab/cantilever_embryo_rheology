
typedef struct {

  Mat A;      

  Vec b;      
  Vec b0;     

  Vec x;      
  Vec x_loc;  

  IS  is_mom;     
  IS  is_cont;    

  IS  is_mom_non_dir;  
  IS  is_cont_non_dir; 

} FullLinearSystem_Stokes_Struct;


PetscErrorCode FullLinearSystem_Stokes_Free(System_Struct *system, FullLinearSystem_Stokes_Struct *stokesFullLS) 
{
  MatDestroy(&stokesFullLS->A);
  VecDestroy(&stokesFullLS->b);
  VecDestroy(&stokesFullLS->x);
  
  ISDestroy(&stokesFullLS->is_mom);
  ISDestroy(&stokesFullLS->is_cont);

  ISDestroy(&stokesFullLS->is_mom_non_dir);
  ISDestroy(&stokesFullLS->is_cont_non_dir);

  if(system->SYSType == IBM)
    VecDestroy(&stokesFullLS->b0);
  return 0; 
}