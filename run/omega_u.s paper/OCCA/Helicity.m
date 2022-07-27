%% Helicity (H) calculation

im1 = @(F) circshift(F, [+1 0]);
jm1 = @(F) circshift(F, [0 +1]);
ip1 = @(F) circshift(F, [-1 0]);
jp1 = @(F) circshift(F, [0 -1]);

A_X = @(F) (F + im1(F)) / 2; %average
D_X = @(F)  F - im1(F); % backwward difference
A_Y = @(F) (F + jm1(F)) / 2;
D_Y = @(F)  F - jm1(F);


dx = g.DXCvec; % dx
dy = g.DYCsc; % dy
A_vor = g.RAZvec; % the area of the vorticity cell

% H on omega_hel
% s_grad_T
s_grad_T_X = A_X(s_hel).*D_X(t_hel); % on u cell
s_grad_T_Y = A_Y(s_hel).*D_Y(t_hel); % on v cell

s_grad_T = (s_grad_T_X + ip1(s_grad_T_X) + s_grad_T_Y + jp1(s_grad_T_Y))./A_vor; % on Vorticity cell

[rs, rt, rsz, rtz] = densjmd95_bsq_second_derivs(s_hel, t_hel, z_hel);

% Thermobaric parameter
tb_hel = (rsz .* rt ./ rs -  rtz); % on tracer cell
tb_hel_vor = (tb_hel + ip1(tb_hel) + jp1(tb_hel) + ip1(jp1(tb_hel)))/4; % on vorticity cell

% Helicity [m^3 kg^-2]
H_hel = s_grad_T.*tb_hel_vor;
