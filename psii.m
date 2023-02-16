function y=psii(u,v)

n=length(v);

y=(2*n)*log(norm(u)^2+norm(v)^2);

for k=1:n
    y=y-log(v(k)^2);
end