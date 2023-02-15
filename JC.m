function Y=JC(N,q,alpha,beta,l,L,R,P)


for a=1:L
    NN(a)=0;
        for p=1:P
            NN(a)=q(p,a)*N(p)+NN(a);
        end
end

for a=1:L
Y(a)=(l(a)^2)/((l(a)*alpha(a)-beta(a)*NN(a))^2+0.01);
end

end