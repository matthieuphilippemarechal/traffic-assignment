function Y=HH(x,q,Delta,Deltaw,alpha,beta,l,L,R,P,DN,Dtilde)


N=x(1:P);
lambd=x(P+1:P+R);
mu=x(P+R+1:2*P+R);

    c=C(N,q,alpha,beta,l,L,P);
% if Dtilde~=-ones(L,1)
%     NN=q'*N;
%     cc=q*(NN./(l'.^2)-Dtilde./l');
%     Y(1:P)=Delta*c+8*cc+Deltaw'*lamb-mu;
% else
    Y(1:P)=Delta*c+Deltaw'*lambd-mu;
% end
   Y(P+1:P+R)=Deltaw*N-DN;

   Y(P+R+1:2*P+R)=N.*mu;

    Y=Y';
end