#!/bin/bash
clear

echo '----------------------------------------------------------------------'
echo ' Compiling on Home PC'

mpicc -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3  -fPIC -Wall -Wwrite-strings -Wno-strict-aliasing -Wno-unknown-pragmas -fstack-protector -fvisibility=hidden -g3    -I/home/mohamad/petsc/petsc_install/include -I/usr/include/valgrind/include  Petsc_IBM_Solver.c  -Wl,-rpath,/home/mohamad/petsc/petsc_install/lib -L/home/mohamad/petsc/petsc_install/lib -Wl,-rpath,/home/mohamad/petsc/petsc_install/lib -L/home/mohamad/petsc/petsc_install/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu/openmpi/lib -L/usr/lib/x86_64-linux-gnu/openmpi/lib -Wl,-rpath,/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9 -Wl,-rpath,/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/lib/x86_64-linux-gnu -L/lib/x86_64-linux-gnu -lpetsc -lcmumps -ldmumps -lsmumps -lzmumps -lmumps_common -lpord -lscalapack -lsuperlu_dist -lflapack -lfblas -lptesmumps -lptscotchparmetis -lptscotch -lptscotcherr -lesmumps -lscotch -lscotcherr -lparmetis -lmetis -lm -lX11 -lstdc++ -ldl -lmpi_usempif08 -lmpi_usempi_ignore_tkr -lmpi_mpifh -lmpi -lgfortran -lm -lgfortran -lm -lgcc_s -lquadmath -lpthread -lrt -lquadmath -lstdc++ -ldl -o Petsc_IBM_Solver
echo ' Compiled Petsc_IBM_Solver'

echo ' Done.'
echo '-----------------------------------------------------------------------'

