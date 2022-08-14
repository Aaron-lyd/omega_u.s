function [z_ntp_im1,z_ntp_ip1,z_ntp_jm1,z_ntp_jp1,T_ntp_im1,T_ntp_ip1,T_ntp_jm1,T_ntp_jp1] = z_T_VENM(S, T, Z, s ,t ,z, OPTS)
% This function is for calculating the NTP depth and the corresponding
% Temperature from the central point (i,j) to its four neighboring casts

% --- Input:
%  S [nk, ni, nj]: practical / Absolute Salinity of the 3D data
%  T [nk, ni, nj]: potential / Conservative Temperature of the 3D data
%  Z [nk, ni, nj] or [nk, 1]: pressure or depth of the 3D data
%  s     [ni, nj]: practical / Absolute Salinity of the surface
%  t     [ni, nj]: potential / Conservative Temperature of the surface
%  z     [ni, nj]: pressure or depth of the surface
%  OPTS  [struct]: options (see "Options" below)
%
%
% --- Output:
%  z_ntp_im1 [ni, nj]: depth of the neutral tangent plane in the x direction minus 1
%  z_ntp_ip1 [ni, nj]: depth of the neutral tangent plane in the x direction plus 1
%  z_ntp_jm1 [ni, nj]: depth of the neutral tangent plane in the y direction minus 1
%  z_ntp_jp1 [ni, nj]: depth of the neutral tangent plane in the y direction plus 1
%  T_ntp_im1 [ni, nj]: temperature on the depth of the neutral tangent plane in the x direction minus 1
%  T_ntp_ip1 [ni, nj]: temperature on the depth of the neutral tangent plane in the x direction plus 1
%  T_ntp_jm1 [ni, nj]: temperature on the depth of the neutral tangent plane in the y direction minus 1
%  T_ntp_jp1 [ni, nj]: temperature on the depth of the neutral tangent plane in the y direction plus 1

% Author(s) : Yandong Lang, Geoff Stanley



% Initialize the functions that will be needed
shift_im1_3D = @(F) circshift(F, [0 +1, 0]); % if z is [lat x lon], this shifts south
shift_ip1_3D = @(F) circshift(F, [0 -1, 0]); % if z is [lat x lon], this shifts south

shift_jm1_3D = @(F) circshift(F, [0 0, +1]);
shift_jp1_3D = @(F) circshift(F, [0 0, -1]);
lead1 = @(x) reshape(x, [1 size(x)]);



[nk,ni,nj] = size(S);
% code generation
ntp_bottle_to_cast_codegen(nk,OPTS);

TOL_DENS = 1e-7;
drdp = 4e-03; % Approximate derivative of in-situ density w.r.t. pressure [kg m^-3 dbar^-1]
X_TOL = TOL_DENS / (drdp * 2); % tolerance in pressure [dbar] during vertical solve.  Factor of 2 for good measure
 
% OPTS.INTERPFN = @ppc_linterp;

SppX = OPTS.INTERPFN(Z, S);
TppX = OPTS.INTERPFN(Z, T);

% ntp-bottle-to-cast in four directions
z_ntp_im1 = nan(size(z));
z_ntp_ip1 = nan(size(z));

z_ntp_jm1 = nan(size(z));
z_ntp_jp1 = nan(size(z));


BotK = squeeze(sum(isfinite(S), 1)); % The number of the non-nan data on the water column
for j = 1:nj
    jm1 = mod(j-2, nj) + 1;
    jp1 = mod(j, nj) + 1;
    for i = 1:ni
        im1 = mod(i-2, ni) + 1; % for periodic date in x direction
        ip1 = mod(i, ni) + 1;
        k = BotK(im1,j);
        z_ntp_im1(i,j) = ntp_bottle_to_cast_mex(SppX(:,:,im1,j), TppX(:,:,im1,j), Z, k, s(i,j), t(i,j), z(i,j), X_TOL);
        k = BotK(ip1,j);
        z_ntp_ip1(i,j) = ntp_bottle_to_cast_mex(SppX(:,:,ip1,j), TppX(:,:,ip1,j), Z, k, s(i,j), t(i,j), z(i,j), X_TOL);
        k = BotK(i,jm1);
        z_ntp_jm1(i,j) = ntp_bottle_to_cast_mex(SppX(:,:,i,jm1), TppX(:,:,i,jm1), Z, k, s(i,j), t(i,j), z(i,j), X_TOL);
        k = BotK(i,jp1);
        z_ntp_jp1(i,j) = ntp_bottle_to_cast_mex(SppX(:,:,i,jp1), TppX(:,:,i,jp1), Z, k, s(i,j), t(i,j), z(i,j), X_TOL);
        
    end
end

% Evaluating the temperature data on the NTP points
T_im1 = shift_im1_3D(T);
T_ip1 = shift_ip1_3D(T);
T_jm1 = shift_jm1_3D(T);
T_jp1 = shift_jp1_3D(T);

T_im1ppX = OPTS.INTERPFN(Z, T_im1);
T_ip1ppX = OPTS.INTERPFN(Z, T_ip1);
T_jm1ppX = OPTS.INTERPFN(Z, T_jm1);
T_jp1ppX = OPTS.INTERPFN(Z, T_jp1);

% T_ntp_im1 = pchipdqn(lead1(z_ntp_im1), Z, T_im1);
% T_ntp_ip1 = pchipdqn(lead1(z_ntp_ip1), Z, T_ip1);
% T_ntp_jm1 = pchipdqn(lead1(z_ntp_jm1), Z, T_jm1);
% T_ntp_jp1 = pchipdqn(lead1(z_ntp_jp1), Z, T_jp1);

% T_venm_x_im1 = interp1qn(lead1(z_venm_x_im1), Z, T_im1);
% T_venm_x_ip1 = interp1qn(lead1(z_venm_x_ip1), Z, T_ip1);
% T_venm_y_jm1 = interp1qn(lead1(z_venm_y_jm1), Z, T_jm1);
% T_venm_y_jp1 = interp1qn(lead1(z_venm_y_jp1), Z, T_jp1);

T_ntp_im1 = ppc_val(Z,T_im1ppX,lead1(z_ntp_im1));
T_ntp_ip1 = ppc_val(Z,T_ip1ppX,lead1(z_ntp_ip1));
T_ntp_jm1 = ppc_val(Z,T_jm1ppX,lead1(z_ntp_jm1));
T_ntp_jp1 = ppc_val(Z,T_jp1ppX,lead1(z_ntp_jp1));


end