% This script compiles the files pardisoReorderLTE.c, pardisoFactorLTE.c
% and pardisoSolveLTE.c on a 64 bit Windows system
% Hints: 
% 1. You first have to install the Intel MKL. This script is written
% for version 9.1.027 and the installation path "C:\Program
% Files\Intel\MKL"
% 2. You have to install a C compiler.

close all;
clear all;

mex -I"C:\Program Files\Intel\MKL\9.1.027\include" -L"C:\Program Files\Intel\MKL\9.1.027\em64t\lib" -lmkl_solver -lmkl_em64t -lmkl_lapack -llibguide40 -v -largeArrayDims pardisoReorderLTE.c
mex -I"C:\Program Files\Intel\MKL\9.1.027\include" -L"C:\Program Files\Intel\MKL\9.1.027\em64t\lib" -lmkl_solver -lmkl_em64t -lmkl_lapack -llibguide40 -v -largeArrayDims pardisoFactorLTE.c
mex -I"C:\Program Files\Intel\MKL\9.1.027\include" -L"C:\Program Files\Intel\MKL\9.1.027\em64t\lib" -lmkl_solver -lmkl_em64t -lmkl_lapack -llibguide40 -v -largeArrayDims pardisoSolveLTE.c


