% call ParDiSo for the solution of a custom sparse system
%
% X = callParDiSo(A, B)
%
% IN: A = system matrix
%     B = right-hand side vectors (column vectors matrix)
%
% OUT: X = solution vectors (column vectors matrix) for each right-hand
%          side vector
%
% Laurent Ntibarikure
function X = callParDiSo(A, B)

solverTime = tic();
fprintf('#> Calling ParDiSo for custom system ...\n');
nbrVectors = size(B,2);

X = zeros(size(A,1),nbrVectors);
for i=1:nbrVectors
  fprintf('##> Solution %d of %d\n', i, nbrVectors);
  
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
  err = pardisoReleaseMemory(Fact.mtype, Fact.iparm, Fact.pt, Fact.A_val, ...
    Fact.A_ia, Fact.A_ja, Fact.ncol);
  if err ~= 0
    error(['Pardiso error during memory release: ' err]);
  end

end

fprintf('#> Solution time %2.4g s.\n', toc(solverTime));