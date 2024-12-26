function [x,error]=pardisoSolveLTE(mtype,iparm,pt,A_val,A_ia,A_ja,ncol, b, ReleaseMemory);
%function [x,error]=pardisoSolveLTE(mtype,iparm,pt,A_val,A_ia,A_ja,ncol, b);
%
% ReleaseMemory is optional and 0 by default
%
% b may contain more than 1 column
%