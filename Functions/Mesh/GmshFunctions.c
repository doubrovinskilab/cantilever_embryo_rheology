#ifndef GMSH_Functions
#define GMSH_Functions


PetscErrorCode GMSHCheckAccuracyOfPartitionedEntitiesDim0(System_Struct *system, PetscInt entities_size, PetscInt *entities_tag, PetscInt *entities_type, 
  PetscInt *entities_group, PetscScalar *x, PetscScalar *y, PetscScalar *z, PetscInt *parent_dim, const PetscInt dim_understudy)
{

  PetscErrorCode ierr = 0;
  PetscInt r, s, i, j;
  PetscScalar *X, *Y, *Z;
  PetscInt S, *C, *c;

  
  MPI_Request Request;
  MPI_Status Status;

  for(r = system->size - 1; r > 0; r--) 
  {

    if ( system->rank == r)
    {
      ierr = MPI_Bcast(&entities_size, 1, MPIU_INT, r, PETSC_COMM_WORLD );      CHKERRMPI(ierr);
      ierr = MPI_Bcast(x,  entities_size, MPIU_SCALAR, r, PETSC_COMM_WORLD );   CHKERRMPI(ierr);
      ierr = MPI_Bcast(y,  entities_size, MPIU_SCALAR, r, PETSC_COMM_WORLD );   CHKERRMPI(ierr);
      ierr = MPI_Bcast(z,  entities_size, MPIU_SCALAR, r, PETSC_COMM_WORLD );   CHKERRMPI(ierr);
    }
    else if (system->rank > r) 
    {
      ierr = MPI_Bcast(&S,            1, MPIU_INT, r, PETSC_COMM_WORLD );       CHKERRMPI(ierr);
      PetscMalloc3(S, &X, S, &Y, S, &Z);
      ierr = MPI_Bcast(X,             S, MPIU_SCALAR, r, PETSC_COMM_WORLD );    CHKERRMPI(ierr);
      ierr = MPI_Bcast(Y,             S, MPIU_SCALAR, r, PETSC_COMM_WORLD );    CHKERRMPI(ierr);
      ierr = MPI_Bcast(Z,             S, MPIU_SCALAR, r, PETSC_COMM_WORLD );    CHKERRMPI(ierr);
      PetscFree3(X, Y, Z);    
    }
    else if (system->rank < r)
    {
      
      ierr = MPI_Bcast(&S,            1, MPIU_INT, r, PETSC_COMM_WORLD );       CHKERRMPI(ierr);
      PetscMalloc4(S, &X, S, &Y, S, &Z, S, &c);
      ierr = MPI_Bcast(X,             S, MPIU_SCALAR, r, PETSC_COMM_WORLD );    CHKERRMPI(ierr);
      ierr = MPI_Bcast(Y,             S, MPIU_SCALAR, r, PETSC_COMM_WORLD );    CHKERRMPI(ierr);
      ierr = MPI_Bcast(Z,             S, MPIU_SCALAR, r, PETSC_COMM_WORLD );    CHKERRMPI(ierr); 

      for (i=0; i<S; i++) 
      { 
        c[i] = 0;
        for (j = 0; j<entities_size; j++) 
          if ((x[j] == X[i]) && (y[j] == Y[i]) && (z[j] == Z[i]) && (parent_dim[j] == dim_understudy))
          { c[i] += 1;
            break; }
      }

      ierr = MPI_Isend( c,  S, MPIU_INT, r, 1000, PETSC_COMM_WORLD, &Request);              CHKERRMPI(ierr);
      ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);

      PetscFree4(X, Y, Z, c);    
    }

    if (system->rank == r)
    {
      PetscMalloc1(entities_size, &C);
      for(i=0; i<entities_size; i++) C[i] = 0;

      
      for (s = 0; s < r; s++) 
      {
        PetscMalloc1(entities_size, &c);
        ierr = MPI_Irecv( c, entities_size, MPIU_INT, s, 1000, PETSC_COMM_WORLD, &Request);   CHKERRMPI(ierr);
        ierr = MPI_Wait(&Request, &Status), MPI_Cancel(&Request), MPI_Request_free(&Request); CHKERRMPI(ierr);
        for(i=0; i<entities_size; i++) C[i] += c[i];
        PetscFree(c);   
      }  
      
      for(i=0; i<entities_size; i++) 
        if (C[i] > 0 && entities_type[i] < 2)
        {
          PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    Correcting the entity '%d' from type '%d' to type '%d' on rank %d\n",entities_tag[i], entities_type[i], 
          entities_type[i]+2, system->rank);
          entities_type[i] += 2;
        }

      PetscFree(C);   
    }
    ierr = PetscSynchronizedFlush(PETSC_COMM_WORLD, PETSC_STDOUT); CHKERRQ(ierr);
    PetscBarrier(NULL);
  }
  
  return 0;
}


PetscErrorCode GMSHCheckVersion(System_Struct *system, const char FLGMSH_File[256])
{ 
  PetscErrorCode ierr;
  char line[256]="";
  char **data;
  PetscInt j;
  PetscScalar gmshVersion;
  PetscBool flg_MeshFormat = PETSC_FALSE;
  FILE *fd;

  fd   = fopen(FLGMSH_File,"r");
  
  if(fd == NULL) 
  { PetscErrorPrintf("GMSH '.msh' file %s is not found.\n", FLGMSH_File);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

  
  while (PETSC_TRUE)
  {
    ReadFileLine(fd, 256,line);

    
    PetscStrncmp(line, "$MeshFormat",           11, &flg_MeshFormat);
    if (flg_MeshFormat)
    {
      ReadFileLine(fd, 256,line);    
      PetscStrToArray(line,' ',&j, &data);
      ierr = PetscOptionsStringToScalar(data[0], &gmshVersion); CHKERRQ(ierr);
      if ( gmshVersion > 4.0) 
        ierr = PetscPrintf(PETSC_COMM_WORLD,"  Gmsh version %.2lf.\n", gmshVersion);
      else{ 
        PetscErrorPrintf("-> INCORRECT Gmsh version %.2lf.\n", gmshVersion);
        SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Incorrect gmsh version. The gmsh version required is > 4.0."); }
      PetscStrToArrayDestroy(j, data);
    }
    else
      continue;
    
    break;
    if (line[0] == 0)  break; 
  }
  
  ierr = fclose(fd);
  if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

  return 0;
}


PetscErrorCode GMSHReadPhysicalNames(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements, const char FLGMSH_File[256])
{ 
  PetscErrorCode ierr;
  char line[256]="";
  char **data;
  char data2[256];
  PetscInt i, j, phys_dim, phys_num;
  PetscBool flg_PhysicalNames = PETSC_FALSE;
  FILE *fd;

  fd   = fopen(FLGMSH_File,"r");
  
  if(fd == NULL) 
  { PetscErrorPrintf("GMSH '.msh' file %s is not found.\n", FLGMSH_File);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

  
  while (PETSC_TRUE)
  {
    ReadFileLine(fd, 256,line);

    
    PetscStrncmp(line, "$PhysicalNames",  14, &flg_PhysicalNames);
    if (flg_PhysicalNames)
    {
      ReadFileLine(fd, 256,line);    
      RemoveNewLine(line,data2);     
      ierr = PetscOptionsStringToInt(data2, &phys_num); CHKERRQ(ierr);

      
      nodes->N_Gr = 0;
      for(i=0; i < phys_num; i++)
      { 
        ReadFileLine(fd, 256,line); 

        PetscStrToArray(line,' ',&j, &data);

        ierr = PetscOptionsStringToInt(data[0], &phys_dim); CHKERRQ(ierr);

        if      (system->Dim == 2 && phys_dim == 1) nodes->N_Gr++;
        else if (system->Dim == 3 && phys_dim == 2) nodes->N_Gr++;

        PetscStrToArrayDestroy(j, data);
      }

      PetscPrintf(PETSC_COMM_WORLD,"  Number of boundary groups: %d\n", nodes->N_Gr);

      elements->N_Gr = nodes->N_Gr;
    }
    else
      continue;
    
    break;
    if (line[0] == 0)  break; 
  }
  
  ierr = fclose(fd);
  if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

  return 0;
}


PetscErrorCode GMSHReadEntities1(System_Struct *system, const char FLGMSH_File[256], PetscInt *entities_size)
{ 
  PetscErrorCode ierr;
  char line[256]="";
  char **data;
  char data2[256];
  PetscInt i, j;
  PetscInt nPar, nGhost;
  PetscBool flg_PartitionedEntities = PETSC_FALSE;
  PetscBool flg_Entities            = PETSC_FALSE;
  FILE *fd;

  
  PetscMPIInt size = system->size;
  

  fd   = fopen(FLGMSH_File,"r");
  
  if(fd == NULL) 
  { PetscErrorPrintf("GMSH '.msh' file %s is not found.\n", FLGMSH_File);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

  
  while (PETSC_TRUE)
  {
    ReadFileLine(fd, 256,line);

    PetscStrncmp(line, "$Entities",              9, &flg_Entities);
    PetscStrncmp(line, "$PartitionedEntities",  20, &flg_PartitionedEntities);

    if ( flg_Entities && size==1)
    {
      ReadFileLine(fd, 256,line);    
      PetscStrToArray(line,' ',&j, &data);

      ierr = PetscOptionsStringToInt(data[0], &entities_size[0]); CHKERRQ(ierr); 
      ierr = PetscOptionsStringToInt(data[1], &entities_size[1]); CHKERRQ(ierr); 
      ierr = PetscOptionsStringToInt(data[2], &entities_size[2]); CHKERRQ(ierr); 
      RemoveNewLine(data[3], data2); 
      ierr = PetscOptionsStringToInt(data2, &entities_size[3]); CHKERRQ(ierr);   
      PetscPrintf(PETSC_COMM_WORLD,"  Number of entities: points=%d, curves=%d, surfaces=%d, volumes=%d\n", 
                            entities_size[0], entities_size[1], entities_size[2], entities_size[3]);
      PetscStrToArrayDestroy(j, data);
    }
    else if  ( flg_PartitionedEntities  && system->size>1)
    {
      ReadFileLine(fd, 256,line);                                    
      PetscStrToArray(line,' ',&j, &data);
      RemoveNewLine(data[0],data2);                                   
      ierr = PetscOptionsStringToInt(data2, &nPar); CHKERRQ(ierr);    
      if( nPar != system->size ){
        PetscErrorPrintf("The GMSH partioned entities is %i different than the number of mpi-cores assigned %i on rank %i.\n", nPar, system->size, system->rank);
        SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! GMSH partitions do not agree with number of mpi-nodes"); }
      PetscStrToArrayDestroy(j, data);

      
      ReadFileLine(fd, 256,line);                                      
      PetscStrToArray(line,' ',&j, &data);  
      RemoveNewLine(data[0],data2);                                     
      ierr = PetscOptionsStringToInt(data2, &nGhost); CHKERRQ(ierr);    
      for(i=0; i< nGhost;i++)   ReadFileLine(fd, 256,line);            
      PetscStrToArrayDestroy(j, data);

      
      ReadFileLine(fd, 256,line);                                      
      PetscStrToArray(line,' ',&j, &data);
      ierr = PetscOptionsStringToInt(data[0],   &entities_size[0]); CHKERRQ(ierr);  
      ierr = PetscOptionsStringToInt(data[1],   &entities_size[1]); CHKERRQ(ierr);  
      ierr = PetscOptionsStringToInt(data[2],   &entities_size[2]); CHKERRQ(ierr);  
      RemoveNewLine(data[3],data2);                                                 
      ierr = PetscOptionsStringToInt(data2,  &entities_size[3]); CHKERRQ(ierr);     
      PetscStrToArrayDestroy(j, data);

      PetscPrintf(PETSC_COMM_WORLD,"  Number of entities:\n");
      PetscSynchronizedPrintf(PETSC_COMM_WORLD,"    rank %d: points=%d, curves=%d, surfaces=%d, volumes=%d\n",system->rank, 
        entities_size[0], entities_size[1], entities_size[2], entities_size[3]); 
      PetscSynchronizedFlush(PETSC_COMM_WORLD,PETSC_STDOUT);
      PetscBarrier(NULL);
    }
    else
      continue;
    
    break;
    if (line[0] == 0)  break; 
  }
  
  ierr = fclose(fd);
  if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

  return 0;
}


PetscErrorCode GMSHReadEntities2(System_Struct *system, const char FLGMSH_File[256], PetscInt *entities_size, 
  PetscInt **entities_tag, PetscInt **entities_type, PetscInt **entities_group, const PetscInt N_Gr)
{ 
  PetscErrorCode ierr;
  char line[PETSC_MAX_PATH_LEN]="";
  char **data;
  PetscInt i,j, k, numPhysicalTags;
  PetscInt numPartition, ownerPartition;
  PetscBool flg_PartitionedEntities = PETSC_FALSE;
  PetscBool flg_Entities            = PETSC_FALSE;
  FILE *fd;
  
  PetscMPIInt size = system->size;
  
  fd   = fopen(FLGMSH_File,"r");
  
  if(fd == NULL) 
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); 
  
  while (PETSC_TRUE)
  {
    ReadFileLine(fd, 256,line);

    PetscStrncmp(line, "$Entities",              9, &flg_Entities);
    PetscStrncmp(line, "$PartitionedEntities",  20, &flg_PartitionedEntities);

    if ( flg_Entities && size==1)
    {
      ReadFileLine(fd, 256,line); 

      
      
      for(j=0; j<4; j++)
        for(i=0; i<entities_size[j]; i++)
        {
          ReadFileLine(fd, 256,line);    
          PetscStrToArray(line,' ',&k, &data);;

          ierr = PetscOptionsStringToInt(data[0], &entities_tag[j][i]); CHKERRQ(ierr); 

          if (j==0)  {ierr = PetscOptionsStringToInt(data[4], &numPhysicalTags); CHKERRQ(ierr);} 
          else       {ierr = PetscOptionsStringToInt(data[7], &numPhysicalTags); CHKERRQ(ierr);} 

          if ( (j==0 && numPhysicalTags == 0) || (j!=0 && numPhysicalTags == 0) )
            entities_group[j][i] = -1; 
          else
            
            if (j==0) {ierr = PetscOptionsStringToInt(data[5], &entities_group[j][i]); CHKERRQ(ierr);}
            else      {ierr = PetscOptionsStringToInt(data[8], &entities_group[j][i]); CHKERRQ(ierr);}
          
          if (j==0 || j==1 || (j==2 && system->Dim == 3)) 
            entities_type[j][i] = 1;
          else 
            entities_type[j][i] = 0;

          PetscStrToArrayDestroy(k, data);
        }

      PetscPrintf(PETSC_COMM_WORLD,"  Entities information is saved\n");
    }
    else if ( flg_PartitionedEntities  && system->size>1)
    {
      for (j=0; j<3; j++) ReadFileLine(fd, 256,line); 

      
      
        PetscInt *parent_dim;
        PetscScalar *x, *y, *z;
        PetscMalloc4(entities_size[0], &x, entities_size[0], &y, entities_size[0], &z, entities_size[0], &parent_dim);
        for(j=0; j<4; j++) 
        {
          for(i=0; i<entities_size[j]; i++)
          {
            ReadFileLine(fd, PETSC_MAX_PATH_LEN,line);    
            PetscStrToArray(line,' ',&k, &data);
            ierr = PetscOptionsStringToInt(data[0], &entities_tag[j][i]); CHKERRQ(ierr);  
            ierr = PetscOptionsStringToInt(data[3], &numPartition); CHKERRQ(ierr);        
            ierr = PetscOptionsStringToInt(data[4], &ownerPartition); CHKERRQ(ierr);      

            if (j==0)  {ierr = PetscOptionsStringToInt(data[numPartition+7],  &numPhysicalTags); CHKERRQ(ierr);} 
            else  { 
              if(k<numPartition+10) printf("Error in %d, line: %s\n",system->rank, line); 
              ierr = PetscOptionsStringToInt(data[numPartition+10], &numPhysicalTags); CHKERRQ(ierr); } 

            if ( (j==0 && numPhysicalTags == 0) || (j!=0 && numPhysicalTags == 0) )
              entities_group[j][i] = -1; 
            else 
              
              if (j==0)  {ierr = PetscOptionsStringToInt(data[numPartition+8],  &entities_group[j][i]); CHKERRQ(ierr);} 
              else       {ierr = PetscOptionsStringToInt(data[numPartition+11], &entities_group[j][i]); CHKERRQ(ierr);} 

            
            if ( entities_group[j][i] > N_Gr )       entities_group[j][i] = -1;

            if (j==0) {
              ierr = PetscOptionsStringToInt(data[1],  &parent_dim[i]);         CHKERRQ(ierr);
              ierr = PetscOptionsStringToScalar(data[numPartition+4],  &x[i]);  CHKERRQ(ierr);
              ierr = PetscOptionsStringToScalar(data[numPartition+5],  &y[i]);  CHKERRQ(ierr);
              ierr = PetscOptionsStringToScalar(data[numPartition+6],  &z[i]);  CHKERRQ(ierr);
            }

            
            if ( ownerPartition == (system->rank+1) ) 
            {
              
              
              if      (j == 0 || j==1)                                    entities_type[j][i] = 1; 
              else if (j == system->Dim-1 && entities_group[j][i] != -1)  entities_type[j][i] = 1; 
              else if (j == system->Dim || j == system->Dim-1)            entities_type[j][i] = 0; 
              else {
                SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_COR,"Error! Unknown entity type for GMSH partition"); }
            }
            else
            {
              
              
              if      (j == 0 || j==1)                                    entities_type[j][i] = 3; 
              else if (j == system->Dim-1 && entities_group[j][i] != -1)  entities_type[j][i] = 3; 
              else if (j == system->Dim || j == system->Dim-1)            entities_type[j][i] = 2; 
              else {
                SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_COR,"Error! Unknown entity type for GMSH partition"); }
            }

            PetscStrToArrayDestroy(k, data);
          }
        }
        PetscFree4(x, y, z, parent_dim);
        PetscBarrier(NULL);
        PetscPrintf(PETSC_COMM_WORLD,"  Partitioned entities information is saved on all ranks\n");
    }
    else
      continue;
    
    break;
    if (line[0] == 0)  break; 
  }
  
  ierr = fclose(fd);
  if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

  return 0;
}


PetscErrorCode GMSHCreateNodes(System_Struct *system, Nodes_Struct *nodes, const char FLGMSH_File[256])
{

  PetscErrorCode ierr;
  FILE *fd;
  char line[256]="";
  char **data;
  char data2[256];
  PetscBool flg_Nodes = PETSC_FALSE;
  PetscInt i,j,k;
  PetscInt LI;                                
  PetscInt *GI;                               
  PetscInt blocks_total, num_nodes_block;     
  PetscInt entity_dim, entity_tag;            

  
  PetscMPIInt size = system->size;
  PetscMPIInt rank = system->rank;      

  PetscPrintf(PETSC_COMM_WORLD,"  Creating nodes from GMSH file\n");

  fd   = fopen(FLGMSH_File,"r");
  
  if(fd == NULL) 
  { PetscErrorPrintf("GMSH '.msh' file %s is not found.\n", FLGMSH_File);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }

  
  while (PETSC_TRUE)
  {
    
    ReadFileLine(fd, 256,line);

    PetscStrncmp(line, "$Nodes", 6, &flg_Nodes);
    if (flg_Nodes)
    {
      ReadFileLine(fd, 256, line);    
      PetscStrToArray(line,' ',&k, &data);
      ierr = PetscOptionsStringToInt(data[0], &blocks_total);         CHKERRQ(ierr);  
      ierr = PetscOptionsStringToInt(data[1], &nodes->N_U_total_loc);  CHKERRQ(ierr); 
      NodesMallocUArraysForAll(nodes, size, rank);
      PetscStrToArrayDestroy(k, data);

      LI = 0;
      
      for(i =0; i < blocks_total; i++ )
      {
        ReadFileLine(fd, 256,line);    
        PetscStrToArray(line,' ',&k, &data);       
        ierr = PetscOptionsStringToInt(data[0], &entity_dim);         CHKERRQ(ierr); 
        ierr = PetscOptionsStringToInt(data[1], &entity_tag);         CHKERRQ(ierr); 
        RemoveNewLine(data[3], data2); 
        ierr = PetscOptionsStringToInt(data2, &num_nodes_block);    CHKERRQ(ierr);  
        PetscMalloc1(num_nodes_block, &GI);  
        
        PetscStrToArrayDestroy(k, data);
        
        for(j =0; j < num_nodes_block; j++ ) 
        { 
          ReadFileLine(fd, 256,line);                                   
          PetscStrToArray(line,' ',&k, &data);
          RemoveNewLine(data[0], data2); 
          ierr = PetscOptionsStringToInt(data2, &GI[j]); CHKERRQ(ierr); 
          GI[j]--;                                                      
          PetscStrToArrayDestroy(k, data);
        } 
        
        for(j =0; j < num_nodes_block; j++ ) 
        {
          ReadFileLine(fd, 256,line);                                    
          PetscStrToArray(line,' ',&k, &data);
          ierr = PetscOptionsStringToScalar(data[0], &nodes->x[LI]);    CHKERRQ(ierr); 
          ierr = PetscOptionsStringToScalar(data[1], &nodes->y[LI]);    CHKERRQ(ierr); 
          RemoveNewLine(data[2], data2); 
          ierr = PetscOptionsStringToScalar(data2, &nodes->z[LI]);      CHKERRQ(ierr);
          nodes->t[LI] = 0;         
          nodes->GI_U[LI] = GI[j];  
          ++LI;
          PetscStrToArrayDestroy(k, data);
        }
        PetscFree(GI);
      }
    }
    else
      continue;
    
    break;
    if (line[0] == 0)  break; 
  }
  
  ierr = fclose(fd);
  if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

  PetscBarrier(NULL);
  PetscPrintf(PETSC_COMM_WORLD,"    Done creating nodes\n");

  return 0;
}


PetscErrorCode GMSHRemoveNodesWithNoElements(System_Struct *system, Elements_Struct *elements, Nodes_Struct *nodes, PetscInt *elem_counter, 
  const PetscInt elem_b_s, const PetscInt elem_i_s, const PetscInt N_U_total_glo, PetscInt *OWNER_LI)
{
  PetscInt i, j, k;
  PetscInt N_del = 0;
  PetscInt *nodes_del;      

  for(i=0; i<nodes->N_U_total_loc; i++) 
    if(elem_counter[i] == 0)
      N_del++;

  if(N_del==0)
  {
    return 0;
  }
  else
  {
    printf("  Rank %d, needs to delete %d nodes\n", system->rank, N_del);
    PetscMalloc1(N_del, &nodes_del); 
    for(i=0; i<N_del;i++) nodes_del[i] = -1;
    
    j=0;
    for(i=0; i<nodes->N_U_total_loc; i++) 
      if(elem_counter[i] == 0){
        nodes_del[j]=i; j++; }
    
    PetscScalar *x_temp, *y_temp, *z_temp;
    PetscInt *t_temp, *GI_U_temp;
    PetscInt sU = nodes->N_U_total_loc;
    PetscMalloc5(sU, &x_temp, sU, &y_temp, sU, &z_temp, sU, &t_temp, sU, &GI_U_temp);
    for(i=0; i<sU; i++)
    { x_temp[i] = nodes->x[i],
      y_temp[i] = nodes->y[i], 
      z_temp[i] = nodes->z[i], 
      t_temp[i] = nodes->t[i],
      GI_U_temp[i] = nodes->GI_U[i];}

    for(i=0; i < elements->N_inter_loc; i++)
      for(j=0; j < elem_i_s; j++)
        elements->inter[i][j] = nodes->GI_U[elements->inter[i][j]];

    for(i=0; i < elements->N_bound_loc; i++)
      for(j=0; j < elem_b_s; j++)
        elements->bound[i][j] = nodes->GI_U[elements->bound[i][j]];

    
    for(PetscInt i=0; i<nodes->N_U_total_loc; i++) PetscFree(nodes->e[i]);
    PetscFree(nodes->e);
    PetscFree6(nodes->x, nodes->y, nodes->z, nodes->t, nodes->GI_U, nodes->LI_UtoP);
    for(PetscInt i=0; i<system->size; i++)  PetscFree(nodes->N_cores[i]);
    PetscFree(nodes->N_cores);

    nodes->N_U_total_loc = nodes->N_U_total_loc - N_del;
    PetscInt sU_new = nodes->N_U_total_loc;
    NodesMallocUArraysForAll(nodes, system->size, system->rank);

    k = 0;
    for(i=0; i<sU_new; i++) 
    { 
      
      for (j=0; j<N_del; j++)
        if (nodes_del[j] == k) k++;

      nodes->x[i]    = x_temp[k];
      nodes->y[i]    = y_temp[k];
      nodes->z[i]    = z_temp[k];
      nodes->t[i]    = t_temp[k];
      nodes->GI_U[i] = GI_U_temp[k];
      nodes->LI_UtoP[i] = -1;
      k++;
    }

    if (k != sU_new+N_del){
      PetscErrorPrintf(" Error not all nodes have been converted.\n");
      SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error in GMSHRemoveNodesWithNoElements"); }
    
    PetscFree(nodes_del);
    PetscFree5(x_temp, y_temp, z_temp, t_temp, GI_U_temp);
    
    for(PetscInt i=0; i<N_U_total_glo; ++i)         OWNER_LI[i] = -1;
    for(PetscInt i=0; i<nodes->N_U_total_loc; ++i)  OWNER_LI[nodes->GI_U[i]] = i;
    
    for(i=0; i < elements->N_inter_loc; i++)
      for(j=0; j < elem_i_s; j++)
        elements->inter[i][j] = OWNER_LI[elements->inter[i][j]];

    for(i=0; i < elements->N_bound_loc; i++)
      for(j=0; j < elem_b_s; j++)
        elements->bound[i][j] = OWNER_LI[elements->bound[i][j]];

  }


  return 0;
}


PetscErrorCode GMSHCreateElements(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements, const char FLGMSH_File[256],
  PetscInt *entities_size, PetscInt **entities_tag, PetscInt **entities_type, PetscInt **entities_group, const PetscInt N_U_total_glo,
  PetscInt *OWNER_LI)
{ 

  PetscErrorCode ierr;
  FILE *fd;
  char line[256]="";
  char **data;
  char data2[256]="";
  PetscBool flg_Elements  = PETSC_FALSE;
  PetscInt i, j, k, b, l, ei, eb, eg, temp;
  PetscInt LI;                                  
  PetscInt GI;                                  
  PetscInt blocks_total, elem_total;            
  PetscInt entity_dim, entity_tag;              
  PetscInt num_elem_block;                      
  PetscInt elemBDim, elemIDim;                  
  PetscInt e;                                   
  PetscInt **elem;                              
  PetscInt elem_t, elem_s, elem_i_s, elem_b_s;  
  PetscInt *elem_type, *elem_tag;
  PetscInt *elem_counter;                       

  const PetscInt Dim = system->Dim;

  const PetscMPIInt size = system->size;
  const PetscMPIInt rank = system->rank;      

  ElementsDefineType(system, elements); 

  const FLUID_ELEMENT_TYPE U_ET = elements->U_elem_type; 

  elem_i_s = elements->U_S_inter;
  elem_b_s = elements->U_S_bound;
  if (U_ET == P1b) elem_i_s = elem_i_s - 1;

  if(Dim == 3) elemBDim = 2, elemIDim = 3;

  PetscMalloc1(nodes->N_U_total_loc, &elem_counter);
  for(i=0; i<nodes->N_U_total_loc; i++) elem_counter[i] = 0;

  PetscPrintf(PETSC_COMM_WORLD,"  Creating elements from GMSH file\n");
  fd   = fopen(FLGMSH_File,"r");
  
  if(fd == NULL) 
  { PetscErrorPrintf("GMSH '.msh' file %s is not found.\n", FLGMSH_File);
    SETERRABORT(PETSC_COMM_WORLD, PETSC_ERR_ARG_WRONG,"Error! File was not found."); }
  
  while (PETSC_TRUE)
  {
    
    ReadFileLine(fd, 256,line);
    PetscStrncmp(line, "$Elements", 9, &flg_Elements);

    if (flg_Elements)
    {
      ReadFileLine(fd, 256, line);    
      
      PetscStrToArray(line,' ',&k, &data);
      ierr = PetscOptionsStringToInt(data[0], &blocks_total); CHKERRQ(ierr); 
      ierr = PetscOptionsStringToInt(data[1], &elem_total);   CHKERRQ(ierr); 
      PetscStrToArrayDestroy(k, data);

      PetscMalloc3(elem_total, &elem_type, elem_total, &elem, elem_total, &elem_tag);  
      for(i = 0; i < elem_total; ++i) PetscMalloc1(elem_i_s, &elem[i]);
      
      elements->N_inter_loc = 0;
      elements->N_bound_loc = 0;
      elements->N_ghost_loc = 0;
            
      e = 0;
      for(i =0; i < blocks_total; i++ )
      {
        ReadFileLine(fd, 256,line);    
        PetscStrToArray(line,' ',&k, &data);

        ierr = PetscOptionsStringToInt(data[0], &entity_dim);         CHKERRQ(ierr); 
        ierr = PetscOptionsStringToInt(data[1], &entity_tag);         CHKERRQ(ierr); 
        RemoveNewLine(data[3], data2); 
        ierr = PetscOptionsStringToInt(data2, &num_elem_block);  CHKERRQ(ierr); 
        PetscStrToArrayDestroy(k, data);

        if ( entity_dim != elemIDim && entity_dim != elemBDim )
        {
          for(j =0; j < num_elem_block; j++ ) {
            ReadFileLine(fd, 256, line);    
            elem_type[e] = -1;                
            e++; }
        }
        else
        {
          if (entity_dim == elemBDim)       
            elements->N_bound_loc = elements->N_bound_loc + num_elem_block, elem_t = 1, elem_s = elem_b_s;
          else if( entity_dim == elemIDim ) 
            elements->N_inter_loc = elements->N_inter_loc + num_elem_block,  elem_t = 0, elem_s = elem_i_s;

          if(elem_t == 1) 
          {
            
            for(b=0; b<entities_size[entity_dim]; b++)
              if(entities_tag[entity_dim][b] == entity_tag)
                break;
            if (entity_tag != entities_tag[entity_dim][b]) 
            { PetscErrorPrintf("Element entity tag %d is not found in rank %d.\n", entity_tag, rank);
              SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Node entity tag not found."); }

            if (entities_type[entity_dim][b]!=1 && size > 1) 
                elem_t = 2, elements->N_bound_loc-=num_elem_block, elements->N_ghost_loc+=num_elem_block;
          }


          for(j =0; j < num_elem_block; j++ ) 
          {
            ReadFileLine(fd, 256,line);    
            PetscStrToArray(line,' ',&k, &data);

            if(j == 0) 
              if( (unsigned) elem_s != k-2 && (unsigned) elem_s != k-1){
                PetscErrorPrintf(" Error size of element is %d while size wanted is %d. Element line: %s\n", (k-2), elem_s, line);
                SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! Wrong element type in InputFile."); }

            for(l=0; l<elem_s; l++) 
            {
              RemoveNewLine(data[l+1], data2); 
              ierr = PetscOptionsStringToInt(data2, &GI); CHKERRQ(ierr); GI--;
              LI = OWNER_LI[GI];
              if(LI == -1){
                PetscErrorPrintf(" Error node with GI=%d was not found in OWNER_LI on rank.\n", GI, rank);
                SETERRABORT(PETSC_COMM_WORLD,PETSC_ERR_ARG_WRONG,"Error! node was not found."); }
              elem[e][l] = LI;        
              elem_counter[LI]++;
            }

            elem_type[e] = elem_t;    

            
            if (elem_t == 2 || elem_t == 0)        elem_tag[e] = -1; 
            else        elem_tag[e] = entities_group[entity_dim][b]; 
            e++;

            PetscStrToArrayDestroy(k, data);
          }
        }
      }
    }
    else
      continue;

    break;
    if (line[0] == 0)  break; 
  }
  
  ierr = fclose(fd);
  if (ierr) SETERRABORT(PETSC_COMM_SELF,PETSC_ERR_SYS,"fclose() failed on file");

  ElementsMallocArrays(system, elements);

  ei = 0, eb = 0, eg = 0;
  for(i=0; i < elem_total; ++i)
  {
    if(elem_type[i] == -1)      
      continue;
    else if (elem_type[i] == 1){ 
      elements->Gr_bound[eb] = elem_tag[i] -1; 
      for (j = 0; j<elem_b_s; ++j)
        {
          elements->bound[eb][j] = elem[i][j];
          
          if (nodes->t[elem[i][j]] == 0 || nodes->t[elem[i][j]] == 2) nodes->t[elem[i][j]]++;
        }
      eb++;}
    else if (elem_type[i] == 0){ 
      for (j = 0; j<elem_i_s; ++j)
        elements->inter[ei][j] = elem[i][j];
      ei++;}
    else if (elem_type[i] == 2) 
      eg++;
  }
  
  if(Dim == 3 && U_ET == P2)
  { 
    for(i=0; i < elements->N_inter_loc; i++ )
    { temp = elements->inter[i][9];
      elements->inter[i][9] = elements->inter[i][8];
      elements->inter[i][8] = temp; }
  }

  for(i=0; i < elements->N_bound_loc; i++) 
    for (j = 0; j<elem_b_s; ++j)
      nodes->t[elements->bound[i][j]] = 1;

  for(i=0; i < nodes->N_U_total_loc; i++ ) 
    if (nodes->t[i] >=2)
      nodes->t[i] -= 2;

  GMSHRemoveNodesWithNoElements(system, elements, nodes, elem_counter, elem_b_s, elem_i_s, N_U_total_glo, OWNER_LI);
  PetscFree(elem_counter);

  if (size == 1) 
  {
    elements->N_inter_glo = elements->N_inter_loc;
    elements->N_bound_glo = elements->N_bound_loc;
  }
  else
  {
    SumIntAcrossRanks(size, rank, &elements->N_inter_glo, elements->N_inter_loc);
    SumIntAcrossRanks(size, rank, &elements->N_bound_glo, elements->N_bound_loc);
  }
  
  for(i = 0; i < elem_total; ++i) PetscFree(elem[i]);
  PetscFree3(elem_type, elem, elem_tag);
  PetscBarrier(NULL);
  PetscPrintf(PETSC_COMM_WORLD,"  Done reading the elements section. Global info: %d internal, %d boundary.\n", elements->N_inter_glo, elements->N_bound_glo);

  return 0;
}


PetscErrorCode GMSHCreateFluidMesh(System_Struct *system, Nodes_Struct *nodes, Elements_Struct *elements)
{
  
  PetscErrorCode ierr;
  PetscLogDouble Clock1, Clock2;
  char FLGMSH_File[256] = "";

  PetscInt N_U_total_glo;  
  PetscInt *OWNER_RANK;    
  PetscInt *OWNER_LI;      
  PetscInt *OWNER_T;       

  
  PetscInt *entities_size, **entities_tag, **entities_type, **entities_group;
  ierr = PetscMalloc4(4, &entities_size, 4, &entities_tag, 4, &entities_type, 4, &entities_group); CHKERRQ(ierr);
  for(PetscInt i=0; i<4; i++) entities_size[i] = 0;


  ierr = PetscBarrier(NULL); CHKERRQ(ierr);
  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"  Computational time: %lf seconds\n\n", (Clock2 - Clock1) );

  
  if(system->size==1) 
    PetscSNPrintf(FLGMSH_File,sizeof(FLGMSH_File),"%s.msh",system->FLMesh_File);
  else                
    PetscSNPrintf(FLGMSH_File,sizeof(FLGMSH_File),"%s_%d.msh",system->FLMesh_File,system->rank+1);

  PetscFixFilename(FLGMSH_File,FLGMSH_File);

  PetscTime(&Clock1);

  PetscPrintf(PETSC_COMM_WORLD,"Reading GMSH '.msh' files on all ranks\n");
  PetscPrintf(PETSC_COMM_WORLD,"-------------------------------------------\n"); 

  GMSHCheckVersion(system, FLGMSH_File);

  GMSHReadPhysicalNames(system, nodes, elements, FLGMSH_File);

  GMSHReadEntities1(system, FLGMSH_File, entities_size);

  for(PetscInt i=0; i<4; i++) 
    PetscMalloc3(entities_size[i], &entities_tag[i], entities_size[i], &entities_type[i], entities_size[i], &entities_group[i]);
  
  GMSHReadEntities2(system, FLGMSH_File, entities_size, entities_tag, entities_type, entities_group, nodes->N_Gr);

  GMSHCreateNodes(system, nodes, FLGMSH_File);

  N_U_total_glo = NodesFindTheMaxGI(system, nodes) + 1;
  PetscMalloc3(N_U_total_glo, &OWNER_RANK, N_U_total_glo, &OWNER_LI, N_U_total_glo, &OWNER_T);
  for(PetscInt i=0; i<N_U_total_glo; i++) OWNER_RANK[i] = -1, OWNER_LI[i] = -1, OWNER_T[i]=-1; 
  for(PetscInt i=0; i<nodes->N_U_total_loc; i++) OWNER_LI[nodes->GI_U[i]] = i;  
  PetscPrintf(PETSC_COMM_WORLD,"    Total global number of nodes is %d\n", N_U_total_glo);

  GMSHCreateElements(system, nodes, elements, FLGMSH_File, entities_size, entities_tag, entities_type, entities_group, N_U_total_glo, OWNER_LI);


  if (N_U_total_glo != NodesFindTheMaxGI(system, nodes) + 1) {
    
    N_U_total_glo = NodesFindTheMaxGI(system, nodes) + 1;
    PetscFree3(OWNER_RANK, OWNER_LI, OWNER_T);
    PetscMalloc3(N_U_total_glo, &OWNER_RANK, N_U_total_glo, &OWNER_LI, N_U_total_glo, &OWNER_T);
    for(PetscInt i=0; i<N_U_total_glo; i++) OWNER_RANK[i] = -1, OWNER_LI[i] = -1, OWNER_T[i]=-1; 
    for(PetscInt i=0; i<nodes->N_U_total_loc; i++) OWNER_LI[nodes->GI_U[i]] = i;  
    PetscPrintf(PETSC_COMM_WORLD,"    Updated total global number of nodes is %d\n", N_U_total_glo);
  }


  MeshFindGhostNodes(system, nodes, elements, N_U_total_glo, OWNER_RANK, OWNER_LI, OWNER_T);

  MeshCreateP1bElements(system, nodes, elements); 

  MeshSaveUBoundaryAndGhostLI(system, nodes, elements, N_U_total_glo, OWNER_RANK, OWNER_LI);

  MeshCreatePNodesFromUNodes(system, nodes, elements);

  MeshReOrganizeNodesGI(system, nodes);

  MeshShareNewGhostGI(system, nodes);

  MeshCheckIfCorrectMesh(system, nodes, elements, N_U_total_glo, OWNER_RANK, OWNER_LI);
  
  PetscFree3(OWNER_RANK, OWNER_LI, OWNER_T);
  for(PetscInt i=0; i<4; i++) {
    ierr = PetscFree3(entities_tag[i], entities_type[i], entities_group[i]); CHKERRQ(ierr);}
  ierr = PetscFree4(entities_size, entities_tag, entities_type, entities_group);CHKERRQ(ierr);

  ierr = PetscBarrier(NULL); CHKERRQ(ierr);
  PetscTime(&Clock2);
  PetscPrintf(PETSC_COMM_WORLD,"  Computational time: %lf seconds\n\n", (Clock2 - Clock1) );

  return 0;
}



#endif