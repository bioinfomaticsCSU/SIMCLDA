function retVal=Qlpl(B,Z,L,qlpl1,qlpl2)
%Qlpl compute approximation of Z_n+1

retVal=qlpl1+L/2*(norm((B-Z),'fro'))^2;
retVal=retVal+trace((B-Z)'*qlpl2);

end