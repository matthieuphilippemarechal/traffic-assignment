clear all

alpha=50*ones(1,19);
%alpha=60*ones(1,19);%+10*round(rand(1,19));
%alpha=[50 60 60 50 50 60 50 50 50 50 60 60 50 60 60 60 50 60 60 ];
% alpha(20)=50;
beta=ones(1,19);
l=10*[0.5,0.1,1,0.3,0.2,0.9,0.2,1,0.5,1,0.2,0.2,0.2,0.2,0.2,0.3,0.3,0.2,0.2];

p(1,:)=[1,7,12,0,0,0];
p(2,:)=[2,4,5,6,12,0];
p(3,:)=[2,4,5,11,15,0];

p(4,:)=[2,4,10,14,15,0];
p(5,:)=[2,9,13,14,15,0];
p(6,:)=[2,9,13,14,18,19];
p(7,:)=[2,9,16,17,19,0];

od(1:7)=1;

p(8,:)=[2,4,5,11,18,0];
p(9,:)=[2,4,10,14,18,0];
p(10,:)=[2,9,13,14,18,0];
p(11,:)=[2,9,16,17,0,0];    

od(8:11)=2;

p(12,:)=[3,4,5,6,12,0];
p(13,:)=[3,4,5,11,15,0];
p(14,:)=[3,4,10,14,15,0];
p(15,:)=[3,9,13,14,15,0];
p(16,:)=[8,13,14,15,0,0];
p(17,:)=[8,13,14,18,19,0];
p(18,:)=[8,16,17,19,0,0];
% p(19,:)=[20,19,0,0,0,0];

od(12:18)=3;

p(19,:)=[3,4,5,11,18,0];
p(20,:)=[3,4,10,14,18,0];
p(21,:)=[3,9,13,14,18,0];
p(22,:)=[8,13,14,18,0,0];
p(23,:)=[8,16,17,0,0,0];
% p(25,:)=[20,0,0,0,0,0];

od(19:23)=4;


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
q2(1:P,1:L)=0;
for j=1:P
    for k=1:M
        if p(j,k)~=0
        Delta(j,p(j,k))=1;
        q(j,p(j,k))=l(p(j,k))/ll(j);
        q2(j,p(j,k))=t(p(j,k))/tt(j);
        end
    end
end

Deltaw(1:R,1:P)=0;

for j=1:P
    Deltaw(od(j),j)=1;
end

Tw=[810;400;450;800];
Tw=[500;500;500;500];
Hcost=[];

[F2,FF2,cF,ccF,lamb2,mu2,err2]=asignacionflujos(Delta,Delta,Deltaw,alpha,beta,l,L,R,P,Tw,od,Hcost);

DN=-Tw.*lamb2;

tol=0.1;
itermax=25;
iter=0;
err(1)=2*tol;
[N,NN,D,cQ,ccQ,Dcrit,TTw,lamb1,mu1,err1(iter+1),nb(iter+1)]=asignacioncantidades(q,Delta,Deltaw,alpha,beta,l,L,R,P,DN);
Dtilde=D;
%q=time_incidence2(cQ,ccQ,Delta,L,P);
while (iter<itermax)&&(err(iter+1)>tol)
qtilde=q;
[q,erri(iter+1),NN]=implicit_time_incidence(Delta,L,P,alpha,beta,l,N);
iter=iter+1;    
qtilde=q;
[N,NN,D,cQ,ccQ,Dcrit,TTw,lamb1,mu1,err1(iter+1),nb(iter+1)]=asignacioncantidades(q,Delta,Deltaw,alpha,beta,l,L,R,P,DN,Dtilde);
FF1=D.*(alpha'-beta'.*D);
F1=N./cQ;
%q=time_incidence2(cQ,ccQ,Delta,L,P);
%err(iter+1)=min([norm(FF1-Delta'*F1),norm(Dtilde-D)]);
err(iter+1)=norm(FF1-Delta'*F1)/norm(FF1);
% k=0;
% while (err(iter+1)>err(iter))&&(k<5)
%     q=0.5*(q+qtilde);
%     [N,NN,D,cQ,ccQ,Dcrit,Tw,lamb1,mu1,err1(iter+1),nb(iter+1)]=asignacioncantidades(qtilde,Delta,Deltaw,alpha,beta,l,L,R,P,DN,Dtilde);
%     FF1=D.*(alpha'-beta'.*D);
%     err(iter+1)=norm(FF1-Delta'*F1)/norm(FF1);
%     k=k+1;
% end
%F1=N./cQ;
%Y=C(N,q,alpha,beta,l,L,P);
%err(iter)=max([norm(Y),err2(iter+1)]);
Dtilde=D;
end

error=norm(Tw-TTw)/norm(Tw);



%Hcost=diseno(D,Dcrit);

% Tw=  [847.3337;
%   549.9993;
%   500.2356;
%   720.6143]



dF=2*norm(FF1-FF2)/(norm(FF1+FF2));
dc=2*norm(ccQ-ccF)/(norm(ccQ+ccF));

promF1=(1/P)*sum(F1);
promF2=(1/P)*sum(F2);

covF1F2=(1/P)*sum((F1-promF1*ones(P,1)).*(F2-promF2*ones(P,1)));
sigmaF1=sqrt((1/P)*sum((F1-promF1*ones(P,1)).^2));
sigmaF2=sqrt((1/P)*sum((F2-promF2*ones(P,1)).^2));

corF1F2=covF1F2/(sigmaF1*sigmaF2);

promcQ=(1/P)*sum(cQ);
promcF=(1/P)*sum(cF);

covcQcF=(1/P)*sum((cQ-promcQ*ones(P,1)).*(cF-promcF*ones(P,1)));
sigmacQ=sqrt((1/P)*sum((cQ-promcQ*ones(P,1)).^2));
sigmacF=sqrt((1/P)*sum((cF-promcF*ones(P,1)).^2));

corcQcF=covcQcF/(sigmacQ*sigmacF);

promFF1=(1/L)*sum(FF1);
promFF2=(1/L)*sum(FF2);

covFF1FF2=(1/L)*sum((FF1-promFF1*ones(L,1)).*(FF2-promFF2*ones(L,1)));
sigmaFF1=sqrt((1/L)*sum((FF1-promFF1*ones(L,1)).^2));
sigmaFF2=sqrt((1/L)*sum((FF2-promFF2*ones(L,1)).^2));

corFF1FF2=covFF1FF2/(sigmaFF1*sigmaFF2);

promccQ=(1/L)*sum(ccQ);
promccF=(1/L)*sum(ccF);

covccQccF=(1/L)*sum((ccQ-promccQ*ones(L,1)).*(ccF-promccF*ones(L,1)));
sigmaccQ=sqrt((1/L)*sum((ccQ-promccQ*ones(L,1)).^2));
sigmaccF=sqrt((1/L)*sum((ccF-promccF*ones(L,1)).^2));

corccQccF=covccQccF/(sigmaccQ*sigmaccF);

if length(Hcost)>0

promFF1crit=1/(length(Hcost))*sum(FF1(Hcost));
promFF2crit=1/(length(Hcost))*sum(FF2(Hcost));
covFF1FF2crit=(1/L)*sum((FF1(Hcost)-promFF1crit*ones(length(Hcost),1)).*(FF2(Hcost)-promFF2crit*ones(length(Hcost),1)));
sigmaFF1crit=sqrt((1/L)*sum((FF1(Hcost)-promFF1crit*ones(length(Hcost),1)).^2));
sigmaFF2crit=sqrt((1/L)*sum((FF2(Hcost)-promFF2crit*ones(length(Hcost),1)).^2));

corFF1FF2crit=covFF1FF2crit/(sigmaFF1crit*sigmaFF2crit);

promccQcrit=1/(length(Hcost))*sum(ccQ(Hcost));
promccFcrit=1/(length(Hcost))*sum(ccF(Hcost));
covccQccFcrit=(1/L)*sum((ccQ(Hcost)-promccQcrit*ones(length(Hcost),1)).*(ccF(Hcost)-promccFcrit*ones(length(Hcost),1)));
sigmaccQcrit=sqrt((1/L)*sum((ccQ(Hcost)-promccQcrit*ones(length(Hcost),1)).^2));
sigmaccFcrit=sqrt((1/L)*sum((ccF(Hcost)-promccFcrit*ones(length(Hcost),1)).^2));

corccQccFcrit=covccQccFcrit/(sigmaccQcrit*sigmaccFcrit);
end
