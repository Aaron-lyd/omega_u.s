function helicity = hel(s, t, z, A_vor)
% calculate the helicity on a surface

% --- Input:
% s [nx,ny]: salinity on the surface
% t [nx,ny]: temperature on the surface
% z [nx,ny]: depth of the surface
% A_vor [nx,ny]: Area if the vorticity cell
%
% --- Output:
% helcity [nx,ny]: helicity on the surface

% define some functions
im1 = @(F) circshift(F, [0 +1 0]);
jm1 = @(F) circshift(F, [0 0 +1]);
ip1 = @(F) circshift(F, [0 -1 0]);
jp1 = @(F) circshift(F, [0 0 -1]);

% avergece and backward difference
A_X = @(F) (F + im1(F)) / 2;
D_X = @(F)  F - im1(F);
A_Y = @(F) (F + jm1(F)) / 2;
D_Y = @(F)  F - jm1(F);

% s_grad_T without divided by the grid distance for the line integral
% around the vorticity cell
s_grad_T_u = A_X(s).*D_X(t); % on u cell
s_grad_T_v = A_Y(s).*D_Y(t); % on v cell

s_grad_T = (-s_grad_T_u + jm1(s_grad_T_u) + s_grad_T_v - im1(s_grad_T_v))./A_vor; % on Vorticity cell

% Thermobaric parameter
[rs, rt, rsz, rtz] = densjmd95_bsq_second_derivs(s, t, z);
tb = (rsz .* rt ./ rs -  rtz); % on tracer cell
tb_vor = (tb + im1(tb) + jm1(tb) + im1(jm1(tb)))/4; % on vorticity cell

helicity = s_grad_T .* tb_vor;

end