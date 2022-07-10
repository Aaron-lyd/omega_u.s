function e_cb = e_cb(K, S, T, Z, s, t, z, t_ip1, t_im1, t_jp1, t_jm1, DXG, DXC, DYG, DYC, Area, integral)
% The capital letters of S, T, Z are 3D data and the lower letters of s, t, z are 2d data

lead1 = @(x) reshape(x, [1 size(x)]);

% calculates grad_n T . grad_n T
if integral
    out = grad_dot_grad(t, t_ip1, t_im1, t_jp1, t_jm1, t, t_ip1, t_im1, t_jp1, t_jm1, DXG, DXC, DYG, DYC, Area);
else
    out = grad_dot_grad_normal(t, t_ip1, t_im1, t_jp1, t_jm1, t, t_ip1, t_im1, t_jp1, t_jm1, DXG, DXC, DYG, DYC, Area);
end 

[rs, rt, ~, ~, rss, rtt, rst, ~] = densjmd95_bsq_second_derivs(s, t, z);
% cabbeling parameter
cb = 2*rst.*rt./rs - rtt - rss.*((rt./rs).^2);

Tppc = ppc_pchip(Z,T); % 3d data on the grid center (g.RC)
Tz = ppc_deriv(Z,Tppc,1); % 3d Tz on the grid center
Tz_2d = ppc_val(Z,Tz,lead1(z)); % interpolate on the 2d surface

Sppc = ppc_pchip(Z,S);
Sz = ppc_deriv(Z,Sppc,1);
Sz_2d = ppc_val(Z,Sz,lead1(z));

% The vertical gradient of the locally referenced potential density
sigma_z = -rs.*Sz_2d - rt.*Tz_2d;

e_cb = K*cb.*out./sigma_z;

end
