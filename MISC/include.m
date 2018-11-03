%%
%%-------------------- whether or not to precondition CG
%%
 CG_prec = 0;
%%
%%-------------------- polynomial degree for chebyshev 
%%
 poldeg = 10; 
%%
%%-------------------- method for diagonalization 
%%  diagmeth == 0 --> Lanczos 1st step and chebyshev filtering thereafter
%%  diagmeth == 1 --> Lanczos all the time
%%  diagmeth == 2 --> Full-Chebyshev subspace iteration first step
%%                    chebyshev filtering thereafter. 
%%
 diagmeth = 0;
%% 
%%-------------------- restart for projection/ mixing 
%%
 restart_dim = 5;
