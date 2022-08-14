function [grad_n_T_x, grad_n_T_y] = grad_n_T_fun(S, T, Z, z, OPTS)



[ni,nj] = size(z);
nk = size(Z,1);

% Just-in-time code generation:
ntp_slope_codegen(nk, ni, nj, isvector(Z));

% Initialize the functions that will be needed
shift_im1_3D = @(F) circshift(F, [0 +1, 0]); % if z is [lat x lon], this shifts south
shift_jm1_3D = @(F) circshift(F, [0 0, +1]);

im1 = @(F) circshift(F, [+1, 0]); % if x is [lat x lon], this shifts south
jm1 = @(F) circshift(F, [0, +1]); % if x is [lat x lon], this shifts west

lead1 = @(x) reshape(x, [1 size(x)]);

interpfn = @ppc_pchip;
SppZ = interpfn(Z, S);
TppZ = interpfn(Z, T);

tolz = OPTS.tolz;
[dzi, dzj] = ntp_slope_mex(SppZ, TppZ, Z, z, tolz, 1, 1);

z_ntp_i = (z + im1(z))/2 - dzi;
z_ntp_j = (z + jm1(z))/2 - dzj;

z_ntp_im1 = (z + im1(z))/2 + dzi;
z_ntp_jm1 = (z + jm1(z))/2 + dzj;


% Evaluating the temperature data on the NTP points
T_im1 = shift_im1_3D(T);
T_jm1 = shift_jm1_3D(T);

T_im1ppX = OPTS.INTERPFN(Z, T_im1);
T_jm1ppX = OPTS.INTERPFN(Z, T_jm1);

T_ntp_i = ppc_val(Z,TppZ,lead1(z_ntp_i));
T_ntp_j = ppc_val(Z,TppZ,lead1(z_ntp_j));

T_ntp_im1 = ppc_val(Z,T_im1ppX,lead1(z_ntp_im1));
T_ntp_jm1 = ppc_val(Z,T_jm1ppX,lead1(z_ntp_jm1));

dx = OPTS.dx;
dy = OPTS.dy;

grad_n_T_x = (T_ntp_i - T_ntp_im1)./dx;
grad_n_T_y = (T_ntp_j - T_ntp_jm1)./dy;
