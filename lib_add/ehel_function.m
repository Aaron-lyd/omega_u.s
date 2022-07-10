function [e_hel, e_hel_rms, product_imh, product_jmh, product_iph, product_jph] = ehel_function(S, T, Z, s, t, z, u, v, DXG, DYG, Area, OPTS)
% Description: Calculate e^hel by using the divergence method
% --- Input:
%  S [nk, ni, nj]: practical / Absolute Salinity of the 3D data
%  T [nk, ni, nj]: potential / Conservative Temperature of the 3D data
%  Z [nk, ni, nj] or [nk, 1]: pressure or depth of the 3D data
%  s     [ni, nj]: practical / Absolute Salinity of the surface
%  t     [ni, nj]: potential / Conservative Temperature of the surface
%  z     [ni, nj]: pressure or depth of the surface
%  u [nk, ni, nj]: horizontal velocity of the 3D data
%  v [nk, ni, nj]: meridional velocity of the 3D data
%  DXG   [ni, nj]: the grid size in the horizontal direction
%  DXG   [ni, nj]: the grid size in the meridional direction
%  Area  [ni, nj]: the area of the grid
%  OPTS  [struct]: options (see "Options" below)
%
%
% --- Output:
%  e_hel     [ni, nj]: e_hel on the omega surface
%  e_hel_rms [ni, nj]: the root mean square of e_hel
% Author(s) : Yandong Lang, Geoff Stanley


ip1 = @(F) circshift(F, [-1, 0]); % if x is [lat x lon], this shifts south
jp1 = @(F) circshift(F, [0, -1]); % if x is [lat x lon], this shifts west
im1 = @(F) circshift(F, [+1, 0]); 
jm1 = @(F) circshift(F, [0, +1]); 

lead1 = @(x) reshape(x, [1 size(x)]);

% z or p on the u and v grid

z_on_U = (z + im1(z)) / 2;
z_on_V = (z + jm1(z)) / 2;

u_on_U = interp1qn(lead1(z_on_U), Z, u);
v_on_V = interp1qn(lead1(z_on_V), Z, v);

% calculating the deoth of the NTP
[z_ntp_im1,z_ntp_ip1,z_ntp_jm1,z_ntp_jp1] = z_T_VENM(S, T, Z, s ,t ,z, OPTS);

% The line integral of each sides of the grid box
product_imh = (z_ntp_im1 - im1(z))/2 .* u_on_U .* DYG./Area;
product_iph = (z_ntp_ip1 - ip1(z))/2 .* ip1(u_on_U .* DYG)./Area;
product_jmh = (z_ntp_jm1 - jm1(z))/2 .* v_on_V .* DXG./Area;
product_jph = (z_ntp_jp1 - jp1(z))/2 .* jp1(v_on_V .* DXG)./Area;

e_hel = (-product_imh + product_iph - product_jmh + product_jph);

e_hel_rms = e_hel - nanmean(e_hel(:));
good = ~isnan(e_hel_rms);
e_hel_rms = rms(e_hel_rms(good));  % RMS of error of omega

end