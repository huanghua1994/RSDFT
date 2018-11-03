function [val]= interp(xs,x,f)
nmax=length(x);
del=x(2)-x(1);
n=fix(xs/del)+1;
f1=f(n);
nm=n+1;
if nm > nmax
    nm=nmax;
end
f2=f(nm);
val=f1+(f2-f1)*(xs-x(n))/del;