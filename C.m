function Y=C(N,q,alpha,beta,l,L,P)

for a=1:L
    NN(a)=0;
        for p=1:P
            NN(a)=q(p,a)*N(p)+NN(a);
        end
end
for a=1:L
Y(a,1)=(l(a)^2)/(l(a)*alpha(a)-beta(a)*NN(a));
end

end