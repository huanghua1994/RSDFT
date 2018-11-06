function [fl, occup] = occupations(ev, T, Nelec, tol)
% function [fl, occup] = occupations(ev, T, Nelec, tol)
% Bisection algorithm for getting Fermi level and occupation factors
% Input:
%   ev    : Eigenvalues computed in diagonalization
%   T     : Temperature
%   Nelec : Total number of electrons
%   tol   : Tolerance of Fermi level error
% Output:
%   fl    : Fermi level
%   occup : Occupation factors
    occup = zeros(size(ev));
    lmax  = ceil(Nelec/2) + 1;
    a = min(ev) - 1;
    b = ev(lmax) + 1; 
    c = (b + a) / 2;

    [ferr_a, ~] = FermiDirac(ev, a, T, Nelec);
    [ferr_b, ~] = FermiDirac(ev, b, T, Nelec);
    [ferr_c, occup] = FermiDirac(ev, c, T, Nelec);

    error = 2 * sum(occup) - Nelec;

    if (ferr_a * ferr_b > tol) 
        c = b;
        occup = ones(1 : lmax)
        disp('[???] In occupations - ferr_a * ferr_b > 0') 
        return
    end 

    iter = 0;
    max_iter = 200;
    while ((abs(error) > tol) && (iter < max_iter))
        iter = iter + 1;
        c = (b + a) / 2;
        [ferr_c, occup] = FermiDirac(ev, c, T, Nelec);
        error = 2 * sum(occup) - Nelec;
        if (ferr_c * ferr_b < 0)
            a = c;
            ferr_a = ferr_c;
        else
            b = c;
            ferr_b = ferr_c;
        end
    end
    
    fl = c;
end