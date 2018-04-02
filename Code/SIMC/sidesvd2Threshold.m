function retVal=sidesvd2Threshold(svdt1,svdt2,svdt3,L,lambda)
%sidesvd2Threshold soft-threshold shrink svd
%function retVal=svdThreshold(U,hatT,B0,L,lambda)

% A=B0-1/L*U'*((U*B0-hatT).*double(hatT>0));
% k=rank(A);
% [L_svd,S_svd,T_svd] = svds(A,k);
% S_svd=diag(max(0,diag(S_svd)-lambda/L));  
% retVal=L_svd*S_svd*T_svd';

A=svdt1-svdt2/L+svdt3/L;
k=rank(A);
[L_svd,S_svd,T_svd] = svds(A,k);
S_svd=diag(max(0,diag(S_svd)-lambda/L));
retVal=L_svd*S_svd*T_svd';

end