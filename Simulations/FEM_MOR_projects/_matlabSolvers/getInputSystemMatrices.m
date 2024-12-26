% Returns the input system matrices of the radiation model
%
% [A, B, X] = getInputSystemMatrices(flag, ports)
%
% IN: flag = 1: for solution, 2: load only
%     nbrPorts = number of ports (WavePorts/LumpedPorts)
% OptIN: critFlag = solves iteratively for each port in order to save
%                   memory
%
% OUT: A = FEM matrix (sparse)
%      B = right-hand side for each port assuming guided fundamental mode
%          (sparse)
%      X = solutions related to each port (full)
%
% Laurent Ntibarikure
function [A, B, X] = getInputSystemMatrices(flag, nbrPorts, critFlag)
% Use ParDiSo (external mex compiled solver) if true or direct internal 
% solver if false
pardisoFlag = true;
% Memory requirements are critical : solves iteratively and save to file
% then load matrices
if nargin>2
  crit = critFlag;
else
  crit = false;
end
%% Input system matrices
if flag == 1 && crit % solves iteratively saving to file then load
  disp('#> Calling ParDiSo ...')
  A = mmread('lte_fileset\systemMatrix.mm');
  fprintf('##> Check available memory for ParDiSo then press a ');
  fprintf('button ...\n    - if not enough press CTRL-C -\n');
  pause
  %%%%% From O.Farle's code for ParDiSo and internal solver
  for i=1:nbrPorts
    fprintf('##> Port %d\n', i);
    b = vectorReader(['lte_fileset\rhs_',num2str(i-1),'.fvec']);
    if pardisoFlag
      % determine matrix type
      if nnz(imag(A))
        Fact.mtype = 6;  % complex symmetric indefinite matrix
      else
        Fact.mtype = -2;  % real symmetric indefinite matrix
      end
      % Fill-reduction analysis and symbolic factorization
      [Fact.iparm Fact.pt Fact.err Fact.A_val Fact.A_ia Fact.A_ja Fact.ncol] = ...
        pardisoReorderLTE(Fact.mtype, A);
      % Numerical factorization
      err = pardisoFactorLTE(Fact.mtype, Fact.iparm, Fact.pt, Fact.A_val, ...
        Fact.A_ia, Fact.A_ja, Fact.ncol);
      if err ~= 0
        error(['Pardiso error during factorization: ' err]);
      end
      % Forward and Backward solve
      [x err] = pardisoSolveLTE(Fact.mtype, Fact.iparm, Fact.pt, ...
        Fact.A_val, Fact.A_ia, Fact.A_ja, Fact.ncol, b, 0);
      if err ~= 0
        error(['Pardiso error during solving: ' err]);
      end
    else
      disp('#> Calling Internal Solver ...')
      [L, U, P, Q] = lu(A);
      x = Q * (U \ (L \ (P * b)));
    end
    if pardisoFlag
      err = pardisoReleaseMemory(Fact.mtype, Fact.iparm, Fact.pt, Fact.A_val, ...
        Fact.A_ia, Fact.A_ja, Fact.ncol);
      if err ~= 0
        error(['Pardiso error during memory release: ' err]);
      end
    end
    clear Fact;
    writeVector(x, ['lte_fileset\sol_',num2str(i-1),'.fvec']);        
  end
  clear b x;
  disp('##> Input System solved, loading Matrices ...');
  B = zeros(size(A,1),nbrPorts);
  X = zeros(size(A,1),nbrPorts);
  for i=1:nbrPorts
    B(:,i) = vectorReader(['lte_fileset\rhs_',num2str(i-1),'.fvec']);
    X(:,i) = vectorReader(['lte_fileset\sol_',num2str(i-1),'.fvec']);
  end
elseif flag == 1 && ~crit % read all the matrices and compute solutions 
  disp('#> Calling ParDiSo ...')
  A = mmread('lte_fileset\systemMatrix.mm');
  B = zeros(size(A,1),nbrPorts);
  X = zeros(size(A,1),nbrPorts);
  fprintf('##> Check available memory for ParDiSo then press a ');
  fprintf('button ...\n    - if not enough press CTRL-C -\n');
  pause
  %%%%% From O.Farle's code for ParDiSo and internal solver
  for i=1:nbrPorts
    fprintf('##> Port %d\n', i);
    B(:,i) = vectorReader(['lte_fileset\rhs_',num2str(i-1),'.fvec']);
    if pardisoFlag
      % determine matrix type
      if nnz(imag(A))
        Fact.mtype = 6;  % complex symmetric indefinite matrix
      else
        Fact.mtype = -2;  % real symmetric indefinite matrix
      end
      % Fill-reduction analysis and symbolic factorization
      [Fact.iparm Fact.pt Fact.err Fact.A_val Fact.A_ia Fact.A_ja Fact.ncol] = ...
        pardisoReorderLTE(Fact.mtype, A);
      % Numerical factorization
      err = pardisoFactorLTE(Fact.mtype, Fact.iparm, Fact.pt, Fact.A_val, ...
        Fact.A_ia, Fact.A_ja, Fact.ncol);
      if err ~= 0
        error(['Pardiso error during factorization: ' err]);
      end
      % Forward and Backward solve
      [X(:,i) err] = pardisoSolveLTE(Fact.mtype, Fact.iparm, Fact.pt, ...
        Fact.A_val, Fact.A_ia, Fact.A_ja, Fact.ncol, B(:,i), 0);
      if err ~= 0
        error(['Pardiso error during solving: ' err]);
      end
    else
      disp('#> Calling Internal Solver ...')
      [L, U, P, Q] = lu(A);
      X(:,i) = Q * (U \ (L \ (P * B(:,i))));
    end

    if pardisoFlag
      err = pardisoReleaseMemory(Fact.mtype, Fact.iparm, Fact.pt, Fact.A_val, ...
        Fact.A_ia, Fact.A_ja, Fact.ncol);
      if err ~= 0
        error(['Pardiso error during memory release: ' err]);
      end
    end
    clear Fact;
    writeVector(X(:,i), ['lte_fileset\sol_',num2str(i-1),'.fvec']);        
  end
  disp('##> Input System solved...')
elseif flag == 2 % load previously solved matrices
  disp('#> Loading Input System matrices ...')
  A = mmread('lte_fileset\systemMatrix.mm');
  B = zeros(size(A,1),nbrPorts);
  X = zeros(size(A,1),nbrPorts);
  for i=1:nbrPorts
    B(:,i) = vectorReader(['lte_fileset\rhs_',num2str(i-1),'.fvec']);
    X(:,i) = vectorReader(['lte_fileset\sol_',num2str(i-1),'.fvec']);
  end
end
A = sparse(A);
B = sparse(B);