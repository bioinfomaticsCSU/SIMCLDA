function [M_recover]=SIMC(M_Omega,Omega_linear,A,B)
%ProSVM implements the algorithm Maxide proposed in [1]
%    Syntax
%       [M_recover,]=Maxide(M_Omega,Omega_linear,A,B,lambda,max_iter)
%    Description
%       Maxide takes,
%           M_Omega       - n\times m, the target matrix with only observed
%                           entries when the unobserved entries are 0 .
%           Omega_linear  - A vector recording the observed positions
%                           in the target matrix. If the (i,j)th entry is 
%                           observed, a value (j-1)*n+i is recorded in Omega_linear.
%           A             - the side information matrix where left singular
%                           vectores lie in. 
%           B             - the side information matrix where right singular 
%                           vectores lie in.
%           lambda        - the regularization parameter
%           max_iter      - maximum number of iterations
%      and returns,
%			M_recover     - the recoverd matrix
%           telapsed      - the training time measures in seconds
%[1] Miao Xu, Rong Jin and Zhi-Hua Zhou. Speed up matrix completion with
%side information: application to multi-label learning. In: NIPS'13.


lambda=1;
max_iter=1000;
M_Omega_linear=full(M_Omega(Omega_linear))';
[n,m]=size(M_Omega);
[row,column]=index2spa(Omega_linear,n);
r_a=size(A,2);
r_b=size(B,2);
L=1;
gamma=2;
Z0=zeros(r_a,r_b);
Z=Z0;
alpha0=1;
alpha=1;
i=0;
convergence=zeros(max_iter,1);
svdt3=A'*M_Omega*B;
AZ0BOmega=xumm(A*Z0,B',row,column);
AZBOmega=AZ0BOmega;
while i<max_iter
    i=i+1;
    Y=Z+alpha*(1/alpha0-1)*(Z-Z0);
    Z0=Z;
    AYBOmega=(1+alpha*(1/alpha0-1))*AZBOmega-(alpha*(1/alpha0-1))*AZ0BOmega;   
    AZ0BOmega=AZBOmega;
    svdt2=A'*(sparse(row,column,AYBOmega,n,m))*B;  
    Z=sidesvd2Threshold(Y,svdt2,svdt3,L,lambda);        
    AZBOmega=xumm(A*Z,B',row,column);   
    qlpl1=norm(AYBOmega-M_Omega_linear,'fro')^2/2;
    qlpl2=svdt2-svdt3;
    DiffL2=norm(AZBOmega-M_Omega_linear,'fro')^2/2;
	
    while DiffL2>Qlpl(Z,Y,L,qlpl1,qlpl2)	% APG
        L=L*gamma;
        Z=sidesvd2Threshold(Y,svdt2,svdt3,L,lambda);        
        AZBOmega=xumm(A*Z,B',row,column);               
        DiffL2=norm(AZBOmega-M_Omega_linear,'fro')^2/2;
    end
    alpha0=alpha;
    alpha=(sqrt(alpha^4+4*alpha^2)-alpha^2)/2;    
    convergence(i,1)=sideobjectCalc(Z,lambda,DiffL2); 
      if i>1
          if abs(convergence(i,1)-convergence(i-1,1))<(1e-5)*convergence(i,1)
                break;
          end
      end
end
M_recover=A*Z*B';


