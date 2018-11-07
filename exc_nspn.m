function [vxc, exc] = exc_nspn(rho, grid_vol) 
% function [vxc, exc] = exc_nspn(rho, grid_vol) 
% LDA, Ceperly-Alder exchange correlation
% Input:
%   rho      : Charge density at given grid points
%   grid_vol : Volume of a gird cell
% Output:
%   vxc : Exchange-correlation values at given grid points
%   exc : Exchange-correlation energy

    g  = -0.2846;  b1 = 1.0529;
    b2 =  0.3334;  c1 = 0.0622;
    c2 =  0.0960;  c3 = 0.0040;
    c4 =  0.0232;  c5 = 0.0192;
    
    one_third = 1/3; 
    a0 = (4 / (9 * pi))^one_third;
    twovpia0 = 2 / (pi * a0);
    p75vpi = 0.75 / pi;
    
    % Ceperly-Alder exchange correlation
    % rs is the local value of the Wigner-Seitz radius
    ngpt = length(rho);
    vxc  = zeros(ngpt, 1);
    exc  = 0;
    for i = 1 : ngpt
        rho_i = rho(i);
        vxc_i = 0;
        if (rho_i > 0) 
            rs = (p75vpi / rho_i)^one_third;
            vxc_i = -twovpia0 / rs;
            exc = exc + 0.75 * vxc_i * rho_i;
            if (rs >= 1)
                sqrs = sqrt(rs);
                ec = g / (1 + b1 * sqrs + b2 * rs);
                vxc_i = vxc_i + ec * ec * (1 + 3.5 * b1 * sqrs * one_third + 4 * b2 * rs * one_third) / g;
            else
                alpha = log(rs);
                ec = c1 * alpha - c2 + (c3 * alpha - c4) * rs;
                vxc_i = vxc_i + ec - (c1 + (c3 * alpha - c5) * rs) * one_third;
            end
            exc = exc + rho_i * ec;
        end
        vxc(i) = vxc_i;
    end 
    
    % Scale the total energy integral by the grid cell volume
    exc  = exc * grid_vol;
end