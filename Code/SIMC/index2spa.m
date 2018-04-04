function [row,column]=index2spa(Omega_linear,n) %
%index2spa compute column and row of the postive sample
%   Usage:  [row,column]=index2spa(Omega_linear,n) 
%	Inputs:
%			Omega_linear: indices of positive samples
%			n:	the number of row
%
%	Outputs:
%			row:	the number of row
%			column: the number of column


if size(Omega_linear,1)==1
    Omega_linear=Omega_linear';
end

row=mod(Omega_linear,n);
row(find(row==0))=n;
column=((Omega_linear-row)/n)';
column=column+1;
row=row';