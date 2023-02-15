clear all

alpha=76.3369*ones(1,7);
beta=0.5663*ones(1,7);
l=[3.350,0.830,0.860,2.230,0.730,0.720,3.450];

p(1,:)=[1,0,0];
p(2,:)=[2,4,5];

od(1:2)=1;

p(3,:)=[3,4,6];
p(4,:)=[7,0,0];

od(3:4)=2;

DN=[193	204]';

[P,M]=size(p);
L=max(max(p));
R=max(od);

t(1:L)=0;
for a=1:L
    t(a)=l(a)/alpha(a);
end



ll(1:P)=0;
tt(1:P)=0;
for j=1:P
for k=1:M
    if p(j,k)~=0
    ll(j)=ll(j)+l(p(j,k));
    tt(j)=tt(j)+t(p(j,k));
    end
end
end

Delta(1:P,1:L)=0;
q(1:P,1:L)=0;
for j=1:P
    for k=1:M
        if p(j,k)~=0
        Delta(j,p(j,k))=1;
        q(j,p(j,k))=t(p(j,k))/tt(j);
        end
    end
end

Deltaw(1:R,1:P)=0;

for j=1:P
    Deltaw(od(j),j)=1;
end


tol=0.0001;
itermax=30;
iter=0;
err(1)=2*tol;
[N,NN,D,cQ,ccQ,Dcrit,Tw,lamb1,mu1,err1(iter+1)]=asignacioncantidades(q,Delta,Deltaw,alpha,beta,l,L,R,P,DN);
Dtilde=D;
q=time_incidence2(cQ,ccQ,Delta,L,P);
while (iter<itermax)&(err(iter+1)>tol)
iter=iter+1;    
qtilde=q;
[N,NN,D,cQ,ccQ,Dcrit,Tw,lamb1,mu1,err1(iter+1)]=asignacioncantidades(qtilde,Delta,Deltaw,alpha,beta,l,L,R,P,DN,Dtilde);
FF1=D.*(alpha'-beta'.*D);
F1=N./cQ;
q=implicit_time_incidence(Delta,L,P,alpha,beta,l,N);
err(iter+1)=min([norm(FF1-Delta'*F1),norm(Dtilde-D)]);
Dtilde=D;
end
FF1=D.*(alpha'-beta'.*D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% q=Delta;
% 
% Hcost=diseno(D,Dcrit);
% 
% % Tw=  [847.3337;
% %   549.9993;
% %   500.2356;
% %   720.6143]
% 
% [F2,FF2,cF,ccF,lamb2,mu2,err2]=asignacionflujos(q,Delta,Deltaw,alpha,beta,l,L,R,P,Tw,od,Hcost);
% 
% DF=densidad(FF2,L,alpha,beta,Hcost)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% b=[1,2;                                              % Estructura de la red
%     1,4;
%     3,4;
%     4,5;
%     5,6;
%     6,7;
%     2,7;
%     3,8;
%     4,8;
%     5,9;
%     6,10;
%     7,11;
%     8,9;
%     9,10;
%     10,11;
%     8,12;
%     12,13;
%     10,13;
%     13,11;
%     3,13];
% 
% s=size(b);
% L=s(1);
% 
% DeltaN(1:max(max(b)),1:L)=0;
% 
% for a=1:L
%     DeltaN(b(a,1),a)=-1;
%     DeltaN(b(a,2),a)=1;
% end


