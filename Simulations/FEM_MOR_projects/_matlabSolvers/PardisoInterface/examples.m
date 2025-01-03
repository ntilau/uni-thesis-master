% This script demonstrates how to use the Pardiso interface with various
% matrix types.

close all;
clear all;

%% Example 1: real symmetric indefinite matrix
A = [7 0 1 0 0 2 7 0;
  0 -4 8 0 2 0 0 0;
  1 8 1 0 0 0 0 5;
  0 0 0 7 0 0 9 0;
  0 2 0 0 5 1 5 0;
  2 0 0 0 1 -1 0 5;
  7 0 0 9 5 0 11 0;
  0 0 5 0 0 5 0 5];
A = sparse(A);
% real symmetric indefinite matrix
mtype = -2; 
% Fill-reduction analysis and symbolic factorization
[iparm pt] = pardisoReorderLTE(mtype, A); 
tic
% Numerical factorization
pardisoFactorLTE(pt, mtype, A, iparm); 
toc
neqns = size(A);
% allocate memory for real solution vector
x = zeros(neqns(1), 1);
b = randn(neqns(1), 1);
% Release all internal memory including factorization
releaseMemory = true;
tic
% Forward and Backward solve
pardisoSolveLTE(pt, mtype, A, iparm, b, x, releaseMemory);
toc
% Test residual
display(norm(A * x - b) / norm(b));


%% Example 2: complex symmetric matrix
dim = 300;
A2 = randn(dim) + 1i * randn(dim);
A2 = A2 + A2.';
A2(5, 5) = 0;
A2(51, 51) = 0;
A2 = sparse(A2);
% complex symmetric matrix
mtype2 = 6;
% Fill-reduction analysis and symbolic factorization
[iparm2 pt2] = pardisoReorderLTE(mtype2, A2);
tic
% Numerical factorization
pardisoFactorLTE(pt2, mtype2, A2, iparm2);
toc
neqns2 = size(A2);
b2 = randn(neqns2(1), 1) + 1i * randn(neqns2(1), 1);
% allocate memory for complex solution vector
x2 = (1 + 1i) * ones(neqns2(1), 1);
% Don't release internal memory: factorization can be reused
releaseMemory = false;
tic
% Forward and Backward solve
pardisoSolveLTE(pt2, mtype2, A2, iparm2, b2, x2, releaseMemory);
toc
% Test residual
display(norm(A2 * x2 - b2) / norm(b2));
% Solve for second right hand side
b2 = randn(neqns2(1), 1) + 1i * randn(neqns2(1), 1);
% allocate memory for complex solution vector
x2 = (1 + 1i) * ones(neqns2(1), 1);
% Release internal memory including factorization
releaseMemory = true;
tic
% Forward and Backward solve
pardisoSolveLTE(pt2, mtype2, A2, iparm2, b2, x2, releaseMemory);
toc
% Test residual
display(norm(A2 * x2 - b2) / norm(b2));



%% Example 3: real and symmetric positive definite matrix
dim = 100;
A3 = rand(dim);
A3 = A3 * A3';
A3 = sparse(A3);
% real and symmetric positive definite matrix
mtype3 = 2;
% Fill-reduction analysis and symbolic factorization
[iparm3 pt3] = pardisoReorderLTE(mtype3, A3);
tic
% Numerical factorization
pardisoFactorLTE(pt3, mtype3, A3, iparm3);
toc
neqns = size(A3);
x3exact = (1 : neqns(1))';
b3 =  A3 * x3exact;
% allocate memory for real solution vector
x3 = zeros(neqns(1),1);
% Release internal memory including factorization
releaseMemory = true;
tic
% Forward and Backward solve
pardisoSolveLTE(pt3, mtype3, A3, iparm3, b3, x3, releaseMemory);
toc
% Test residual
display(norm(A3 * x3 - b3) / norm(b3));
% Test accuracy of solution vector in infinity norm
display(norm(x3exact - x3, inf) / norm(x3exact, inf));

