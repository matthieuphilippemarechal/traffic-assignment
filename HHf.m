function Y=HHf(x,Delta,Deltaw,alpha,beta,l,L,R,P,Tw)


F=x(1:P);
lamb=x(P+1:P+R);
mu=x(P+R+1:2*P+R);

    c=Cfapprox(F,Delta,alpha,beta,l,L,P);

    Y(1:P)=Delta*c+Deltaw'*lamb-mu;

   Y(P+1:P+R)=Deltaw*F-Tw;

   Y(P+R+1:2*P+R)=F.*mu;

    Y=Y';
end