 function [W, ritzv] = chefsi1(Vin, ritzv, deg, nev, H) 
%   function [W, ritzv] = chefsi(Vin, ritzv, deg, nev, H) 
%   Vin   -- basis vector
%           
%   ritzv -- vector containing previous retz values
%   deg   -- polynomial degree
%   nev   -- number of occupied states
%   H     -- The hamiltonian matrix 
% out:
%   ritzv -- contains the new ritz values
%   W     -- the filtered vectors
%------------------------------------------------------------------
 [n,n2] = size(Vin);
 G = zeros(n2,n2);
%% 
 if (nev >= n2) 
     error('Vin should have more than nev columns') 
 end
%%
%% call Lanczos with nev = 6 and m=6 to get the upper interval bound
%%  
  [Dum1, Dum2, upperb] = lan(H, 6, randn(n,1), 6, 0.0) ; 
%%
%% use max of previous eigenvalues for lower bound
%% 
 lowerb  = max(ritzv) ; 
%% 
%% use the smallest eigenvalue estimate for scaling.
%% 
  lam1    = min(ritzv);
%--------------------
  if (lowerb > upperb), error('bounds are wrong'); end 
%
  W = ch_filter(H, Vin, deg, lam1, lowerb, upperb); 
%
% orthonormalize the basis, should be replaced with better method
% for real computations
%
  [W, R]=qr(W,0);   
%
%-------------------- Rayleigh-ritz projection - 
%                      compute Hhat = W'*H*W     
  Vin = H*W;
  for j=1:n2
      for i=1:j
          G(i,j) = Vin(:,i)'*W(:,j);
          G(j,i) = G(i,j);
      end
  end 
%-------------------- diagonalize
  [Q,D] = eig(G);
  W  = W*Q;
  ritzv = diag(D);
%-------------------- finished 

