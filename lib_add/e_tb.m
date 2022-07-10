function e_tb = e_tb(K, S, T, Z, s, t, z, t_ip1, t_im1, t_jp1, t_jm1, z_ip1, z_im1, z_jp1, z_jm1, DXG, DXC, DYG, DYC, Area, integral)
% The capital letters of S, T, Z are 3D data and the lower letters of s, t, z are 2d data

lead1 = @(x) reshape(x, [1 size(x)]);

% calculates grad_n T . grad_n P

if integral
    out = grad_dot_grad(t, t_ip1, t_im1, t_jp1, t_jm1, z, z_ip1, z_im1, z_jp1, z_jm1, DXG, DXC, DYG, DYC, Area);
else
    out = grad_dot_grad_normal(t, t_ip1, t_im1, t_jp1, t_jm1, z, z_ip1, z_im1, z_jp1, z_jm1, DXG, DXC, DYG, DYC, Area);
end

[rs, rt, rsz, rtz, ~, ~, ~, ~] = densjmd95_bsq_second_derivs(s, t, z);

% Thermobaric parameter, since the pressure uses dbar as the unit, we
% multiply 1e-4 to make it to be unit of Pa.
tb = (rsz .* rt ./ rs -  rtz);


Tppc = ppc_pchip(Z,T); % 3d data on the grid center (g.RC)
Tz = ppc_deriv(Z,Tppc,1); % 3d Tz on the grid center
Tz_2d = ppc_val(Z,Tz,lead1(z)); % interpolate on the 2d surface

Sppc = ppc_pchip(Z,S);
Sz = ppc_deriv(Z,Sppc,1);
Sz_2d = ppc_val(Z,Sz,lead1(z));

% The vertical gradient of the locally referenced potential density
sigma_z = -rs.*Sz_2d - rt.*Tz_2d;

e_tb = K*tb./sigma_z.*out;
end