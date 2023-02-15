function [q,err,NN]=implicit_time_incidence(Delta,L,P,alpha,beta,l,N)

itermax=2000;
tol=0.01;
err=2*tol;
iter=1;
NN=(alpha./(2*beta)).*l;
NN=NN';
while (iter<itermax)&&(err>tol)
    GG=G(Delta,L,P,alpha,beta,l,N,NN);
    JGG=JG(Delta,L,P,alpha,beta,l,N,NN);
    d=JGG\(-GG);
    k=0;
    kmax=10;
    eps=0.8;
    GGG=G(Delta,L,P,alpha,beta,l,N,NN+eps^k*d);
    while ((min(NN+eps^k*d)<0)||(max(NN+eps^k*d-((alpha./beta).*l)')>0)||(norm(GGG)^2>(1-2*eps^(k+1))*norm(GG)^2))&&(k<kmax)
        k=k+1;
        GGG=G(Delta,L,P,alpha,beta,l,N,NN+eps^k*d);        
    end
    NN=NN+eps^k*d;
    iter=iter+1;
    err=norm(GGG);
end
tt=(l.^2)./(l.*alpha-beta.*NN');

t(1:P)=0;
for p=1:P
   t(p)=sum(Delta(p,:).*tt); 
end

for p=1:P
    for a=1:L
        q(p,a)=(tt(a)/t(p))*Delta(p,a);
    end
end
%err=max([norm(GG),-min(NN+eps^k*d),max(NN+eps^k*d-((alpha./beta).*l)')]);

end