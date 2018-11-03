  function [y] = ch_filter(A, x, deg, lam1, low, high) 
% function [y] = ch_filter(A, x, deg, lam1, low, high) 
%  --> Apply chebyshev filter to x. 
%  A    = matrix
%  x    = vector (s) to be filtered 
%  lam1 = estimate of lowest eigenvalue - for scaling
%         purposes only [rough estimate OK] 
%  [low, high] = interval to be damped. 
% 
  e = (high - low)/2;
  c = (high+low)/2;
  sigma1 = e/(lam1 - c);  
  sigma  = sigma1;
%%-------------------- degree 1 term
    y = (A*x - c*x) .* (sigma1/e);
%%-------------------- loop to degree
    for i = 2: deg
       sigma_new = 1. /(2./sigma1 - sigma);
       t1 = 2*sigma_new/e; t2 = sigma*sigma_new;
       ynew = (A*y - c*y)*t1 - t2*x;
       x = y;             
       y = ynew;          
       sigma = sigma_new;  
    end
    
