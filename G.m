function y=G(Delta,L,P,alpha,beta,l,N,NN)

tt=(l.^2)./(l.*alpha-beta.*NN');

t(1:P)=0;
for p=1:P
   t(p)=sum(Delta(p,:).*tt); 
end

y(1:L,1)=0;
for a=1:L
   y(a,1)=NN(a)-tt(a)*sum(Delta(:,a).*N./t'); 
end

end