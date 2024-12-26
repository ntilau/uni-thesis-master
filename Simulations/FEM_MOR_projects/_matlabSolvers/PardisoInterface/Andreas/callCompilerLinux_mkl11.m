% This script compiles the files pardisoReorderLTE.c, pardisoFactorLTE.c
% and pardisoSolveLTE.c on a 64 bit Windows system
% Hints: 
% 1. You first have to install the Intel MKL. This script is written
% for version 9.1.027 and the installation path "/opt/intel/mkl"
% 2. This script was tested with the gcc compiler.

%close all;
%clear all;
%%
intelpath='/opt/intel/Compiler/11.0/081';
mklpath=[intelpath '/mkl'];

%libs='-lmkl_solver -lmkl_lapack -lmkl_em64t -lmkl -lvml -lguide -lpthread'; % fuer mkl 9

%'mkl_vml_mc'
mkllibs={'mkl_gf_lp64' 'mkl_gnu_thread' 'mkl_core' 'mkl_solver_lp64' 'mkl_em64t'};

libs=[' -L"' mklpath '/lib/em64t/" -L"' intelpath '/lib/intel64/" '];
if(0) %dynamic linking  
  for s=mkllibs
    libs=[libs ' -l' s{:}];
  end  
else %static linking
  for s=mkllibs
    libs=[libs ' ' mklpath '/lib/em64t/lib' s{:} '.a'];
  end
end

libs=[libs ' -lguide -lgomp -lpthread'];
%libs=[libs ' -liomp5 -lpthread'];
%
%%
eval(['mex -I"' mklpath '/include" ' libs ' -v -largeArrayDims pardisoReorderLTE.c']);
eval(['mex -I"' mklpath '/include" ' libs ' -v -largeArrayDims pardisoFactorLTE.c']);
eval(['mex -I"' mklpath '/include" ' libs ' -v -largeArrayDims pardisoSolveLTE.c']);
eval(['mex -I"' mklpath '/include" ' libs ' -v -largeArrayDims pardisoReleaseMemory.c']);

%testen ob's alles korrekt gelinkt wurde
try
  pardisoReorderLTE
catch
  disp(lasterr)
end
try
  pardisoSolveLTE
catch
  disp(lasterr)
end
try
  pardisoFactorLTE
catch
  disp(lasterr)
end
try
  pardisoReleaseMemory
catch
  disp(lasterr)
end