function Y=JHHf(x,q,Delta,Deltaw,alpha,beta,l,L,R,P,Hcost)



F=x(1:P);
lamb=x(P+1:P+R);
mu=x(P+R+1:2*P+R);

jc=JCfapprox(F,q,alpha,beta,l,L,R,P);

Y(1:2*P+R,1:2*P+R)=0;

Y(1:P,1:P)=Delta*diag(jc)*q';
       Y(1:P,P+1:P+R)=Deltaw';
       Y(1:P,P+R+1:2*P+R)=-eye(P);

Y(P+1:P+R,1:P)=Deltaw;
Y(P+R+1:2*P+R,1:P)=diag(mu);
Y(P+R+1:2*P+R,P+R+1:2*P+R)=diag(F);
    
    


end