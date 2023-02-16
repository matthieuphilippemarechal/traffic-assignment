function y=Cf(q,F,P,L,alpha,beta,l,Hcost)



FF=q'*F;

k=1;
for i=1:L
    if k<=length(Hcost)
   if i==Hcost(k)
       d(i)=Dens(FF(i),alpha(i),beta(i),2);
       y(i)=l(i)/(alpha(i)-beta(i)*d(i)+0.01);
       k=k+1;
   else
       d(i)=Dens(FF(i),alpha(i),beta(i),1);
       y(i)=l(i)/(alpha(i)-beta(i)*d(i)); 
   end
    else
            d(i)=Dens(FF(i),alpha(i),beta(i),1);
       y(i)=l(i)/(alpha(i)-beta(i)*d(i));    
    end
   if FF(i)>(alpha(i)^2)/(4*beta(i))
    y(i)=y(i)+100*(FF(i)-(alpha(i)^2)/(4*beta(i)));
    end
end

y=y';
