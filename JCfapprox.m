function Y=JCfapprox(F,Delta,alpha,beta,l,L,R,P)


for a=1:L
    FF(a)=0;
        for p=1:P
            FF(a)=Delta(p,a)*F(p)+FF(a);
        end
end

for a=1:L
    if alpha(a)==50
Y(a)=0.0054674137*0.0005330200*exp(0.0054674137*FF(a));
    else
        Y(a)=0.004422205*l(a)*0.000309373*exp(0.004422205*FF(a));
    end
end

end