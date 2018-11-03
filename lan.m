 function [W, lam, bound] = lan(B, nev, v, m, tol)
%% function [W, lam] = lan(B, nev, v, m)
%% B   = matrix;
%% nev = number of wanted egenvalues;
%% v   = initial vector
%% m   = number of steps
%% tol = tolerance for stopping
%% can do with (full) reorthogonalization (reorth=1)
%% or no -- reorthogonalization (reorth=0) 
 reorth = 0;
%%
n = size(B,1);
v = v/norm(v); 
v1 = v;
v0 = zeros(n,1);
k = 0 ; 
bet = 0; 
ll = 0;
%--------------------pre-allocating
VV   = zeros(n,m);
Tmat = zeros(m+1,m+1);
tr   = zeros(1,m);
%%-------------------- main lanczos loop 
while (k < m)
   k   = k+1;
   VV(:,k) = v1;
   v   =  B*v1 - bet*v0;
   alp = v1'*v;
   v   = v-alp*v1 ;
%%-----------------------------
%% reorth  -- test for reorth. needed! 
%%-----------------------------
   if (reorth) 
     t = VV(:,1:k)'*v ;
     v = v - VV(:,1:k)*t ;
   end
%%-------------------- normalize and store v1
   bet = norm(v); 
   v0 = v1;
   v1 = v/bet;
%%-------------------- update tmat
   Tmat(k,k)   = alp;
   Tmat(k+1,k) = bet;
   Tmat(k,k+1) = bet;
   NTest  = min(5*nev,m) ;     %% when to   start testing
%%-------------------- tr, ll,  == for plotting
   indx = [];
if ( ( (k >= NTest) && (mod(k,10) == 0 )) || k == m )
      rr  = eig(Tmat(1:k,1:k));
%%
%%-------------------- upper bound used by chef_si -- not used otherwise. 
%%    bound = max(abs(rr)) + bet ; 
      bound = max(abs(rr)) + bet;   %% *max(abs(X1(:,k)));
      [rr, indx]  = sort(rr) ;      %% sort increasingly
      tr1 = sum(rr(1:nev));
      ll = ll+1;
      tr(ll) = tr1;
end
%%
%% stopping criterion based on sum of eigenvalues.
%% make sure this is is well-converged!
%%
   if (ll>1 && (abs(tr(ll)-tr(ll-1)) < tol*tr(ll-1)))
      break
   end
%% end - big while loop
 end
%%-------------------- save e.values and compute e.vectors 
   [X1,rr]  = eig(Tmat(1:k,1:k));
   rr = diag(rr);
   [rr, indx]  = sort(rr) ;      %% sort increasingly
%%
   if (reorth) 
      lam = rr(1:mev);
      Y = X1(:,indx(1:nev));
      W = VV(:,1:k)*Y;
   else
%%-------------------- if no reorth. to get meaningful
%%-------------------- eigenvectors.. get twice as many
%%-------------------- eigenvectors and use RayleighRitz       
    mev = min(nev*3,k);
    lam = rr(1:mev);
    Y = X1(:,indx(1:mev));
    W = VV(:,1:k)*Y;
%% add an orthogonalization step 
   [W,G] = qr(W,0) ;
   Vin = B*W;
   for j=1:mev
       for i=1:j
           G(i,j) = Vin(:,i)'*W(:,j);
           G(j,i) = G(i,j);
      end
   end 
   [X1,rr]  = eig(G);
   rr = diag(rr);
   [rr, indx]  = sort(rr) ;      %% sort increasingly
   lam = rr(1:nev);
   Y = X1(:,indx(1:nev));
   W = W*Y;
 end   
%%----------------------------------------------------------- 
