clear all

if(0)
  load test512
  s0=0;
  A=sys.K+s0*sys.D+s0^2*sys.M;
else
  load K1023;
  A=K+K';
end

%b = eye(size(A,1), 3);
p=amd(A);
tic;R=chol(A(p,p));toc



%%
m=2;
A=spdiags((1:3)',0,3,3)+1i*speye(3);
b=ones(size(A,1), m)+0i*eye(size(A,1), m);

%
% mtype  2 real and symmetric positive definite
%       -2 real and symmetric indefinite
%        6 complex and symmetric
%
% Error Information
%   0   No error.
%  -1   Input inconsistent.
%  -2   Not enough memory.
%  -3   Reordering problem.
%  -4   Zero pivot, numerical fact. or iterative refinement problem.
%  -5   Unclassified (internal) error.
%  -6   Preordering failed (matrix types 11, 13 only).
%  -7   Diagonal matrix problem.
%  -8   32-bit integer overflow problem

mtype = 6; 


% Fill-reduction analysis and symbolic factorization
[iparm,pt,err,A_val,A_ia,A_ja,ncol] = pardisoReorderLTE(mtype, A); 
tic
% Numerical factorization
err=pardisoFactorLTE(mtype,iparm,pt,A_val,A_ia,A_ja,ncol);
fprintf(' error = %i\n',err); 
toc
neqns = size(A);
% allocate memory for real solution vector
% x = zeros(size(b));
tic
% Forward and Backward solve
releasememory=0;
[x,err]=pardisoSolveLTE(mtype,iparm,pt,A_val,A_ia,A_ja,ncol, b, releasememory);
fprintf(' error = %i\n',err);
toc
% Test residual
display(norm(A * x - b) / norm(b));

error=pardisoReleaseMemory(mtype,iparm,pt,A_val,A_ia,A_ja,ncol);
