function y=JCf(F,Delta,P,L,alpha,beta,l,Hcost)

for a=1:L
    FF(a)=0;
        for p=1:P
            FF(a)=Delta(p,a)*F(p)+FF(a);
        end
end

y(1:L)=0;


k=1;
for i=1:L
    if k<=length(Hcost)
   if i==Hcost(k)
       d(i)=Dens(FF(i),alpha(i),beta(i),2);
       y(i)=(beta(i)*l(i))/((alpha(i)-2*beta(i)*d(i))*((alpha(i)-beta(i)*d(i)+0.01)^2));
       k=k+1;
   else
       d(i)=Dens(FF(i),alpha(i),beta(i),1);
       y(i)=(beta(i)*l(i))/((alpha(i)-2*beta(i)*d(i))*((alpha(i)-beta(i)*d(i))^2));
   end
    else
        d(i)=Dens(FF(i),alpha(i),beta(i),1);
       y(i)=(beta(i)*l(i))/((alpha(i)-2*beta(i)*d(i))*((alpha(i)-beta(i)*d(i))^2));
    end
end  