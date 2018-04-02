function retValue=sideobjectCalc(B,lambda,DiffL2)
%sideobjectCalc obejective function
s=svd(B);
retValue=lambda*sum(s);
%DiffL2=(norm((U*B-hatT).*double(hatT>0),'fro'))^2/2
retValue=retValue+DiffL2;

end