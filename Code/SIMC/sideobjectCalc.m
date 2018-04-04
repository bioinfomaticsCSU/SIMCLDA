function retValue=sideobjectCalc(B,lambda,DiffL2)
%sideobjectCalc calculate the minimum value of the obejective function

s=svd(B);
retValue=lambda*sum(s);
retValue=retValue+DiffL2;

end