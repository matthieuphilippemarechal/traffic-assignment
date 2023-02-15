function y=gradpsi(u,v)

n=length(u);
m=length(v);
y=((2*m)/log(norm(u)^2+norm(v)^2))*[u;v];
y(n+1:n+m,1)=y(n+1:n+m,1)-2*(1./v);
end