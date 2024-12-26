% This script compiles the files pardisoReorderLTE.c, pardisoFactorLTE.c
% and pardisoSolveLTE.c on a 64 bit Windows system
% Hints: 
% 1. You first have to install the Intel MKL. This script is written
% for version 9.1.027 and the installation path "/opt/intel/mkl"
% 2. This script was tested with the gcc compiler.

close all;
clear all;

mex -I"/opt/intel/mkl/9.1.023/include" -L"/opt/intel/mkl/9.1.023/lib/em64t/" -lmkl_solver -lmkl_lapack -lmkl_em64t -lmkl -lvml -lguide -lpthread  -v -largeArrayDims pardisoReorderLTE.c
mex -I"/opt/intel/mkl/9.1.023/include" -L"/opt/intel/mkl/9.1.023/lib/em64t/" -lmkl_solver -lmkl_lapack -lmkl_em64t -lmkl -lvml -lguide -lpthread  -v -largeArrayDims pardisoFactorLTE.c
mex -I"/opt/intel/mkl/9.1.023/include" -L"/opt/intel/mkl/9.1.023/lib/em64t/" -lmkl_solver -lmkl_lapack -lmkl_em64t -lmkl -lvml -lguide -lpthread  -v -largeArrayDims pardisoSolveLTE.c
