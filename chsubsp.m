 function [W, ritzv] = shsubsp(deg, nev, H) 
%   function [W, ritzv] = shsubsp(deg, nev, H) ;
%
%   deg   -- polynomial degree
%   nev   -- number of occupied states
%   H     -- The hamiltonian matrix 
% out:
%   ritzv -- contains the new ritz values
%   W     -- the approximate invariant subspace
%------------------------------------------------------------------
  Lanc_steps   = max(3*nev,100);
%%%changed Max_out to 20 from 10
%%-------------------- ask for a very moderate precision 
%%                     This is a slow method! 
  Energ_tol    = 0.025;
  Max_out_iter = 10;
  n = size(H,1) ;
%%
%% call Lanczos with Lanc_steps steps for getting the upper interval bound
%%  
  [W, ritzv, upperb] = lan(H, nev, randn(n,1), Lanc_steps, 0.0) ; 
%%
%%-------------------- use previous eigenvalues for lower bound
%% 
  tr0     = sum(ritzv);
%%
%%--------------------outer iterastion 
%%
for it = 1:Max_out_iter 
%%-------------------- use the smallest eigv. estimate for scaling.
   lam1    = min(ritzv);
%--------------------
   lowerb  = max(ritzv) ; 
   if (lowerb > upperb), error('bounds are wrong'); end 
%% 
   W = ch_filter(H, W, deg, lam1, lowerb, upperb); 
%
%-------------------- Rayleigh-ritz projection. 
% orthonormalize the basis, should be replaced with better method
% for real computations
%
   [W, R] = qr(W,0);   
%-------------------- compute Hhat = W'*H*W     
   Vin = H*W;
   [n,n2] = size(Vin);
   for j=1:n2
      for i=1:j
          G(i,j) = Vin(:,i)'*W(:,j);
          G(j,i) = G(i,j);
      end
   end 
%%-------------------- ritz values/vectors 
   [Q, D] = eig(G);
   ritzv  = diag(D);
   [ritzv indx]  = sort(ritzv) ; %%% sort increasingly
   W      = W*Q(:,indx) ;
   tr1    = sum(ritzv(1:nev)) ;
   if (abs(tr1 - tr0) < Energ_tol*abs(tr1))
     break;
   end
   tr0 = tr1;
end
 
%-------------------- 

