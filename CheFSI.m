function [psi_out, rv_out] = CheFSI(H, psi, m, rv)
% function [psi_out, rv_out] = CheFSI(H, psi, m, rv)
% Chebyshev-filtered subspace iteration for SCF 
% Input:
%   H   : Hamiltonian matrix
%   psi : Current wave function
%   m   : Degree of Chebyshev polynomial (optional, default is 8)
%   rv  : Previous Ritz values
% Output:
%   psi_out : Filtered new wave function
%   rv_out  : New Ritz values
    if (nargin == 3), m = 8; end

    % 1. Compute the bounds of current H's eigenvalues
    [lambda, f_2norm] = Lanczos_KSteps(H, 6);
    a0    = min(rv);
    b_low = max(rv);
    b_up  = max(lambda) + f_2norm;
    
    % 2. Perform Chebyshev filtering on previous basis to get new basis
    psi_out = Chebyshev_Filter(H, psi, m, a0, b_low, b_up);
    
    % 3. Orthonormalize new basis
    psi_out = orth(psi_out);
    
    % 4. Perform the Rayleigh-Ritz step
    % (1) Compute Hbar
    Hbar = psi_out' * H * psi_out;
    % (2) Compute H*V = V*diag(D), D has H's non-increasingly ordered Ritz values
    [V, D] = eig(Hbar);
    D = diag(D);
    % Don't know why we cannot use the following 3 statements in RSDFT
    %[~, index] = sort(D, 'descend');
    %V = V(:, index);
    %D = D(index);
    % (3) Rotate the basis
    psi_out = psi_out * V;  
    rv_out = D;
end