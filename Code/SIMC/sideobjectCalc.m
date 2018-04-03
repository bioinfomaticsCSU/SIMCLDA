function retValue=sideobjectCalc(B,lambda,DiffL2)
%sideobjectCalc calculate the minimum value of the obejective function
%   Usage:  retValue=sideobjectCalc(B,lambda,DiffL2)
%	Inputs:
%			B: 		convex approximation solution of Z
%			lambda:	the regularization parameter used to balance the
%			nuclear norm and approximation error
%			DiffL2: approximation error
%
%	Outputs:
%			retValue: the minimum value of the obejective function

s=svd(B);
retValue=lambda*sum(s);
%DiffL2=(norm((U*B-hatT).*double(hatT>0),'fro'))^2/2
retValue=retValue+DiffL2;

end