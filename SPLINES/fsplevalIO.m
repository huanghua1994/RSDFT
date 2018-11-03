  function [y j_out] = fsplevalIO(z,c,d,xi,x,j_in) 
%% [y j_out] = fsplevalB(z,c,d,xi,x,j_in) 
%% evaluates free spline function at x. 
%% x a scalar (only). 
%% z, c, d, as output from fspline.. 
%% j_in  = input value for interval to try first
%%         this interval is [xi(jin) xi(j_in+1)]
%% j_out =  output interval for next call
%%
%% algorithm tests if x is in [x_[j_in] x_[j_in+1] ]
%% if yes -- then computes the value -- if not
%% then finds the correct interval.
%%------------------------------------------ 
 n = length(xi); 
%%-------------------- input j_in is incorrect-- restart:
if (j_in <1 || j_in > n-1) 
     j_in = 1;
end
%%-------------------- x is outside interval -- bring it to
%%                     one of the boundary points
if (x < xi(1)) 
     x = xi(1);
     j_in = 1;
   elseif (x > xi(n))
     x = xi(n);
     j_in = n-1;
end
j_out = j_in;
%%-------------------- if point not in input interval 
%%                     do Binary search
if (not(xi(j_in)<=x && x<=xi(j_in+1)))
    [j_out fflag]=binary_search(xi,x);
    if (fflag ~= 0) 
        error(' SPLINE ERROR [ in binary search ] ')
    end 
end

t1  = xi(j_out+1) - x; 
t2  = x - xi(j_out);
h_j_out = xi(j_out+1) - xi(j_out);
y   = t1*(z(j_out)*t1*t1/h_j_out+c(j_out))+...
      t2*(z(j_out+1)*t2*t2/h_j_out+d(j_out));
