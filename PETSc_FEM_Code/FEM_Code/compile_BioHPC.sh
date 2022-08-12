#!/bin/bash
clear
module purge
module add shared slurm openmpi/gcc/64/3.1.1 
export PATH=$PATH:/work/biophysics/s203893/petsc_2/valgrind-3.17.0/bin
export PATH=$PATH:/work/biophysics/s203893/petsc_2/petsc_install/bin

echo '----------------------------------------------------------------------'
echo ' Compiling on BioHPC'

mpicc -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3  -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3    -I/endosome/work/biophysics/s203893/petsc_2/petsc/include -I/endosome/work/biophysics/s203893/petsc_2/petsc/arch-linux2-c-debug/include -I/work/biophysics/s203893/petsc_2/petsc_install/include -I/work/biophysics/s203893/petsc_2/valgrind-3.17.0/include  Petsc_IBM_Solver.c  -Wl,-rpath,/endosome/work/biophysics/s203893/petsc_2/petsc/arch-linux2-c-debug/lib -L/endosome/work/biophysics/s203893/petsc_2/petsc/arch-linux2-c-debug/lib -Wl,-rpath,/work/biophysics/s203893/petsc_2/petsc_install/lib -L/work/biophysics/s203893/petsc_2/petsc_install/lib -Wl,-rpath,/cm/shared/apps/openmpi/gcc/64/3.1.1-20200326/lib -L/cm/shared/apps/openmpi/gcc/64/3.1.1-20200326/lib -Wl,-rpath,/cm/shared/apps/slurm/16.05.8/lib64 -L/cm/shared/apps/slurm/16.05.8/lib64 -Wl,-rpath,/cm/shared/apps/gcc/5.4.0/lib/gcc/x86_64-unknown-linux-gnu/5.4.0 -L/cm/shared/apps/gcc/5.4.0/lib/gcc/x86_64-unknown-linux-gnu/5.4.0 -Wl,-rpath,/cm/shared/apps/gcc/5.4.0/lib64 -L/cm/shared/apps/gcc/5.4.0/lib64 -Wl,-rpath,/cm/shared/apps/slurm/16.05.8/lib64/slurm -L/cm/shared/apps/slurm/16.05.8/lib64/slurm -Wl,-rpath,/cm/shared/apps/gcc/5.4.0/lib -L/cm/shared/apps/gcc/5.4.0/lib -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu_dist -lflapack -lfblas -lptesmumps -lptscotchparmetis -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lparmetis -lmetis -ltriangle -lm -lX11 -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lrt -lquadmath -lstdc++ -ldl -o Petsc_IBM_Solver
echo ' Compiled Petsc_IBM_Solver'

echo ' Done.'
echo '-----------------------------------------------------------------------'

