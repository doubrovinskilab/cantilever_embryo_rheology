#ifndef Misc_H
#define Misc_H

#include <errno.h>

PetscErrorCode SumScalarAcrossRanks(const PetscMPIInt size, const PetscMPIInt rank, PetscScalar *sum, const PetscScalar add)
{
  /* Does a scalar global sum across all rank */
  PetscErrorCode ierr;
  *sum = 0;
  PetscScalar recv = 0;
  for (PetscInt r=0; r< size; r++)
  {
    if (r == rank) // send add
      recv = add;

    ierr = MPI_Bcast(&recv, 1, MPIU_SCALAR, r, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 

    *sum += recv;
  }

  return(0);
}


PetscErrorCode SumIntAcrossRanks(const PetscMPIInt size, const PetscMPIInt rank, PetscInt *sum, const PetscInt add)
{
  /* Does an int global sum across all rank */
  PetscErrorCode ierr;
  *sum = 0;
  PetscInt recv = 0;
  for (PetscInt r=0; r< size; r++)
  {
    if (r == rank) // send add
      recv = add;

    ierr = MPI_Bcast(&recv, 1, MPIU_INT, r, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 

    *sum += recv;
  }

  return(0);
}


PetscErrorCode MaxIntAcrossRanks(const PetscMPIInt size, const PetscMPIInt rank, PetscInt *global_max, const PetscInt local_max)
{
  /* finds a max int across all rank */
  PetscErrorCode ierr;
  *global_max = -10000000;
  PetscInt recv = 0;
  for (PetscInt r=0; r< size; r++)
  {
    if (r == rank) // send add
      recv = local_max;
    ierr = MPI_Bcast(&recv, 1, MPIU_INT, r, PETSC_COMM_WORLD ); CHKERRMPI(ierr); 
    if (recv > *global_max)
      *global_max = recv;
  }

  return(0);
}


void OMPGetLoopLimits(const PetscInt N, const PetscInt OMP_N_threads, const PetscInt OMP_Rank, 
  PetscInt *i_start, PetscInt *i_end)
{
  /* Divides a loop */
  PetscInt di = (PetscInt) N/OMP_N_threads;
  *i_start = OMP_Rank*di;
  *i_end = (OMP_Rank+1)*di;
  if (OMP_N_threads-1 == OMP_Rank)
    *i_end = N;
}



PetscErrorCode ReadFileLine(FILE *fp,size_t len,char string[])
{
  char *ptr = fgets(string, len, fp);
  if (!ptr) {
    string[0] = 0;
    if (!feof(fp)) {
      SETERRQ1(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Error reading from file: %d", errno);
      SETERRABORT(PETSC_COMM_SELF, PETSC_ERR_FILE_READ, "Error reading from file!");}
  }
  return(0);
}

PetscErrorCode RemoveNewLine(char string_in[], char string_out[])
{
  PetscInt i;
  char **data;
  PetscStrToArray(string_in,'\n',&i, &data);    
  PetscStrncpy(string_out, data[0], 256);
  PetscStrToArrayDestroy(i, data);
  PetscStrToArray(string_out,'\r',&i, &data);   
  PetscStrncpy(string_out, data[0], 256);
  PetscStrToArrayDestroy(i, data);
  return(0);
}

PetscErrorCode RemoveCharAtEnd(char string_in[], char string_out[], char string_remove)
{
  PetscInt i;
  char **data;
  PetscStrToArray(string_in, string_remove,&i, &data);    
  PetscStrncpy(string_out, data[0], 256);
  PetscStrToArrayDestroy(i, data);
  return(0);
}

PetscErrorCode RemoveCharAtBegin(char string_in[], char string_out[], char string_remove)
{
  PetscInt i;
  char **data;
  PetscStrToArray(string_in,string_remove,&i, &data);    
  PetscStrncpy(string_out, data[1], 256);
  PetscStrToArrayDestroy(i, data);
  return(0);
}

PetscErrorCode RemoveCharAtEndAndNewLine(char string_in[], char string_out[], char string_remove)
{
  PetscInt i;
  char **data;
  PetscStrToArray(string_in,string_remove,&i, &data);    
  PetscStrncpy(string_out, data[0], 256);
  PetscStrToArrayDestroy(i, data);
  PetscStrToArray(string_out,'\n',&i, &data);            
  PetscStrncpy(string_out, data[0], 256);
  PetscStrToArrayDestroy(i, data);
  PetscStrToArray(string_out,'\r',&i, &data);            
  PetscStrncpy(string_out, data[0], 256);
  PetscStrToArrayDestroy(i, data);
  return(0);
}



PetscErrorCode RemoveExtraSpaces(char* string) {
  PetscInt i,j;
  size_t lent;
  PetscStrlen(string, &lent);
  PetscInt len = (PetscInt) lent;
  for(i=0; i<len; i++)
  {
    if(string[0]==' '){
      for(i=0; i<(len-1); i++)
      string[i] = string[i+1];
      string[i] = '\0';
      len--;
      i = -1;
      continue;
    }
    if(string[i]==' ' && string[i+1]==' '){
      for(j=i; j<(len-1); j++){
         string[j] = string[j+1];
      }
      string[j] = '\0';
      len--;
      i--;
    }
  }
  return(0);
}



void swap(PetscInt* a, PetscInt* b)
{
  PetscInt t = *a;
  *a = *b;
  *b = t;
}
 
PetscInt partition (PetscInt arr[], PetscInt low, PetscInt high)
{
  PetscInt pivot = arr[high]; 
  PetscInt i = (low - 1); 

  for (PetscInt j = low; j <= high - 1; j++)
  {
      
      if (arr[j] < pivot)
      {
          i++; 
          swap(&arr[i], &arr[j]);
      }
  }
  swap(&arr[i + 1], &arr[high]);
  return (i + 1);
}
 

void quickSort(PetscInt arr[], PetscInt low, PetscInt high)
{
  if (low < high)
  {
      PetscInt pi = partition(arr, low, high);
      quickSort(arr, low, pi - 1);
      quickSort(arr, pi + 1, high);
  }
}

#endif