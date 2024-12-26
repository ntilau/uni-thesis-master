load matrix.mat
A = sparse(i,j,s);
mtype = -2 ;
[iparm, pt] = pardisoinit(mtype);
iparm (3) = 1;
% perform reordering only in the first step
pardisoreorder(pt, mtype, A, iparm);
%  numerical factorization
tic
   pardisofactor(pt, mtype, A, iparm);
toc
neqns = size(A);
x = ones(neqns(1),1);
b =  A*x;
x = zeros(neqns(1),1);
%  solve
tic
   pardisosolve(pt, mtype, A, iparm, b, x);
toc
norm(A*x-b)/norm(b)


