function max_ev = Lanczos_MaxEigVal(H, k)
% Estimate the upper bound of matrix H's eigenvalues by k-step Lanczos
% k is optional, default is 6
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
    max_ev = norm(T, 2) + norm(f, 2);
end