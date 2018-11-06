function [psi, rv] = SCF_Step1_CheFSI(H, n_occ, m)
% function [psi, rv] = SCF_Step1_CheFSI(H, n_occ, m)
% First SCF step diagonalization using CheFSI
% Input:
%   H     : Hamiltonian matrix
%   n_occ : Number of occupied orbitals
%   m     : Degree of Chebyshev polynomial (optional, default is 8)
% Output:
%   psi   : Initial approximate invariant subspace of wave functions
%   rv    : Initial Ritz values
    if (nargin == 2), m = 8; end
    
    % 1. Compute the bounds of current H's eigenvalues
    [lambda, f_2norm] = Lanczos_KSteps(H, 6);
    beta  = 0.5;
    a0    = min(lambda);
    b_up  = max(lambda) + f_2norm;
    b_low = beta * min(lambda) + (1 - beta) * max(lambda);
    
    % 2. Set psi as a random matrix with s >= n_occ column vectors
    n   = size(H, 1);
    s   = n_occ + 10;
    psi = rand(n, s);
    
    % 3. Perform CheFSI several times to get the initial psi
    maxiter = 4;
    for iter = 1 : maxiter
        % (1) Apply Chebyshev filter to psi
        psi = Chebyshev_Filter(H, psi, m, a0, b_low, b_up);
        % (2) Orthonormalize psi
        psi = orth(psi);
        % (3) Compute Rayleigh quotient
        G = psi' * H * psi;
        % (4) Compute eigen-decomposition of G
        [V, D] = eig(G);
        % (5) Sort the Ritz values in increasing order
        [R, index] = sort(diag(D));
        % (6) Reset b_low and a0
        a0    = min(R);
        b_low = max(R);
    end
    % Perform subspace rotation and save the Ritz values
    psi = psi * V(:, index);
    rv  = R;
end