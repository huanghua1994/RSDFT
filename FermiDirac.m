function [ferr, occup] = FermiDirac(ev, mu, T, Nelec) 
% function [ferr, occup] = FermiDirac(ev, mu, T, Nelec) 
% Input:
%   ev    : Eigenvalues computed in diagonalization
%   mu    : Chemical potential
%   T     : Temperature
%   Nelec : Total number of electrons
% Output:
%   ferr  : Error of Fermi level
%   occup : Occupation factors
    kB    = 6.33327186e-06;  % Boltzmann¡¯s constant
    kBT   = kB * T;
    t     = 1 + exp((ev - mu) ./ kBT) ; 
    spin  = 1;
    occup = spin ./ t;
    ferr  = sum(occup) - Nelec/2; 
end