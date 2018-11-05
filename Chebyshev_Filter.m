function Y = Chebyshev_Filter(H, X, m, a0, a, b)
% Filter vectors by a m-degree Chebyshev polynomial on a given interval
% Input:
%   H    : Hamiltonian matrix
%   X    : Vectors to be filtered
%   m    : Degree of Chebyshev polynomial
%   a0   : Estimated lowest eigenvalue, for scaling
%   a, b : Interval [a, b] is to be damped
% Output: 
%   Y    : Filtered vectors
    e = (b - a) / 2;
    c = (b + a) / 2;
    sigma = e / (a0 - c);
    gamma = 2 / sigma;
    Y = (sigma / e) .* (H * X - c .* X);
    for i = 2 : m
        sigma1 = 1 / (gamma - sigma);
        Y1 = (2 * sigma1 / e) .* (H * Y - c .* Y) - (sigma * sigma1) .* X;
        X = Y;
        Y = Y1;
        sigma = sigma1;
    end
end