function [F,FF,cF,ccF,lamb,mu,err]=asignacionflujos(Delta,Deltaw,alpha,beta,l,L,R,P,Tw)


iter=0;
itermax=4500;
tol=10^(-3);
err=2*tol;
x=5*ones(2*P+R,1);
bbeta=0.5;
JH=JHHf(x,Delta,Deltaw,alpha,beta,l,L,R,P);


while (iter<itermax)&&(err>tol)
    H=HHf(x,Delta,Deltaw,alpha,beta,l,L,R,P,Tw);
    JH=JHHf(x,Delta,Deltaw,alpha,beta,l,L,R,P);
rank(JH);
d=(JH)\(-H+[0*ones(P+R,1);((0.5/P)*ones(1,P)*H(P+R+1:2*P+R))*ones(P,1)]);
jj=0;
mF=min(x(1:P)+d(1:P));
mmu=min(x(P+R+1:2*P+R)+d(P+R+1:2*P+R));
m=min([mF,mmu]);
while (m<=-0.1)&&(jj<10)
    jj=jj+1;

mF=min(x(1:P)+(bbeta^jj)*d(1:P));
mmu=min(x(P+R+1:2*P+R)+(bbeta^jj)*d(P+R+1:2*P+R));
m=min([mF,mmu]);
end

if m>-0.1
HP(:,1)=HHf(x+(bbeta^jj)*d,Delta,Deltaw,alpha,beta,l,L,R,P,Tw);

m=-psii(HP(1:P+R),HP(P+R+1:2*P+R))+psii(H(1:P+R),H(P+R+1:2*P+R))+0.2*(bbeta^jj)*gradpsi(H(1:P+R),H(P+R+1:2*P+R))'*JH(1:2*P+R,:)*d;
while (m<-0.1)&&(jj<5)
    jj=jj+1;
    HP(:,1)=HHf(x+(bbeta^jj)*d,Delta,Deltaw,alpha,beta,l,L,R,P,Tw);
    m=-psii(HP(1:P+R),HP(P+R+1:P+2*R))+psii(H(1:P+R),H(P+R+1:P+2*R))+0.2*(bbeta^jj)*gradpsi(H(1:P+R),H(P+R+1:2*P+R))'*JH(1:2*P+R,:)*d;
end
end
x=x+bbeta^jj*d;
iter=iter+1;
err=norm(H);
H=HHf(x,Delta,Deltaw,alpha,beta,l,L,R,P,Tw);
e(iter)=norm(H)+norm(min([x(1:P,1),0*ones(P,1)]))+norm(min([x(P+R+1:2*P+R,1),0*ones(P,1)]));
% x(1:P)=max(x(1:P),0*ones(P,1));
% x(P+R+1:2*P+R)=max(x(P+R+1:2*P+R),0*ones(P,1));
end
F=x(1:P);
lamb=x(P+1:P+R);
mu=x(P+R+1:2*P+R);
ccF=Cfapprox(F,Delta,alpha,beta,l,L,P);
cF=Delta*ccF;
FF=Delta'*F;








end