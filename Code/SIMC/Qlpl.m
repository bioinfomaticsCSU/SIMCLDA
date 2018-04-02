function retVal=Qlpl(B,Z,L,qlpl1,qlpl2)
%qlpl2=U'*((U*Z-hatT).*double(hatT>0));

retVal=qlpl1+L/2*(norm((B-Z),'fro'))^2;
%temp=U'*((U*Z-hatT).*double(hatT>0));
retVal=retVal+trace((B-Z)'*qlpl2);

end