function y=JG(Delta,L,P,alpha,beta,l,N,NN)

tt=(l.^2)./(l.*alpha-beta.*NN');
dtt=(beta.*l.^2)./((l.*alpha-beta.*NN').^2);

t(1:P)=0;
for p=1:P
   t(p)=sum(Delta(p,:).*tt); 
end

y(1:L,1:L)=0;
for a=1:L
   y(a,a)=1-dtt(a)*sum(Delta(:,a)'.*N'.*(t-tt(a)*ones(1,P))./(t.^2)); 
   for aa=1:L
       if a~=aa
           y(a,aa)=dtt(aa)*tt(a)*sum(Delta(:,a)'.*Delta(:,aa)'.*N'./(t.^2));
       end
   end
end

end