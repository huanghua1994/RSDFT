function [x, its] = pcg (A, rhs, x0, m, tol, PRE, precfun)
%% solves A x = rhs by pcg 
%%------------------------------------------ 
n = size(A,1); 
x = x0; 
r = rhs - A * x;
if (nargin >5) 
  z = feval(precfun,PRE,r);
else 
  z = r;
end
p = z ;
ro1 = z' * r; 
tol1 = tol*tol*ro1; 
%%
its = 0 ;
while (its < m && ro1 > tol1) 
  its = its+1; 
  ro = ro1; 
  ap = A * p; 
  alp = ro / ( ap'* p ) ;
  x = x + alp * p ;
  r = r - alp * ap; 
  if (nargin >5) 
    %%
    %% Unpreconditioned case 
    %%
    z = feval(precfun,PRE,r);
  else 
    z = r; 
  end
  %%
  ro1= z'*r;
  bet = ro1 / ro ;
  p = z + bet * p; 
end

 
