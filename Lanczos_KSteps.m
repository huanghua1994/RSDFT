function [lambda, f_2norm] = Lanczos_KSteps(H, k)
% Estimating max / min eigenvalue bounds using k-step Lanczos
% Input:
%   H : Target matrix
%   k : Number of Lanczos steps, (optional, default is 6)
% Output:
%   lambda  : Eigenvalues of tridiagonal matrix T
%   f_2norm : 2-norm of vector f
    if (nargin == 1), k = 6; end;
    
    T = zeros(k);
    n = size(H, 1);
    v = rand(n, 1);
    v = v ./ norm(v, 2);
    f = H * v;
    a = f' * v;
    f = f - a .* v;
    T(1, 1) = a;
    for j = 2 : k
        v0 = v;
        b  = norm(f, 2);
        v  = f ./ b;
        f  = H * v;
        f  = f - b .* v0;
        a  = f' * v;
        f  = f - a .* v;
        T(j, j - 1) = b;
        T(j - 1, j) = b;
        T(j, j) = a;
    end
    lambda  = eig(T);
    f_2norm = norm(f, 2);
end