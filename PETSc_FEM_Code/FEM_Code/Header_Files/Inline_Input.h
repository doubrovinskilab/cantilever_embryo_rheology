
ierr = PetscOptionsGetBool(NULL,NULL,"-print_nodes",&print_nodes,NULL);                                   CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_elements",&print_elements,NULL);                             CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_stokesVariables", &print_stokesVariables, NULL);             CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_solid", &print_solid, NULL);                                 CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_mat_A",&print_mat_A,NULL);                                   CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_vec_b",&print_vec_b,NULL);                                   CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_mat_Ar",     &print_mat_Ar,NULL);                            CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_mat_Pr",     &print_mat_Pr,NULL);                            CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_vec_br",     &print_vec_br,NULL);                            CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_vec_b_dir",  &print_vec_b_dir,NULL);                         CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-check_symmetry_Ar",&check_symmetry_Ar,NULL);                       CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_dir_V_U", &print_dir_V_U,NULL);                              CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_dir_V_P", &print_dir_V_P,NULL);                              CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_non_dir_V_U", &print_non_dir_V_U,NULL);                      CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-print_non_dir_V_P", &print_non_dir_V_P,NULL);                      CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-write_solid_fe", &write_solid_fe,NULL);                      CHKERRQ(ierr);
ierr = PetscOptionsGetBool(NULL,NULL,"-write_fluid_mesh", &write_fluid_mesh,NULL);                        CHKERRQ(ierr);

ierr = PetscOptionsGetString(NULL,NULL,"-in",system.SYSInputFile,sizeof(system.SYSInputFile),&flg);               CHKERRQ(ierr);
if (!flg) SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_USER_INPUT, "Must indicate an ASCII input file with the -in option.");

write_fluid_N_dtime = 0;
write_fluid_info = PETSC_FALSE;
ierr = PetscOptionsGetInt(NULL,NULL,"-write_fluid_info", &write_fluid_N_dtime, &flg); CHKERRQ(ierr);
if (flg && write_fluid_N_dtime > 0) 
  write_fluid_info = PETSC_TRUE;
else  
  write_fluid_info = PETSC_FALSE, write_fluid_N_dtime = 0;
