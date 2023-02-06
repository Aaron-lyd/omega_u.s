function [s, Df, s_rms, Df_rms, e_hel, e_hel_rms, d] = s_ehel_df(z, S, T, Z, U, V, wrap, OPT)

% functions
im1 = @(F) circshift(F, [+1, 0]); % if x is [lat x lon], this shifts south
jm1 = @(F) circshift(F, [0, +1]); % if x is [lat x lon], this shifts west

ip1 = @(D) circshift(D, [-1 0]);
jp1 = @(D) circshift(D, [0 -1]);

A_X_f = @(F) (F + ip1(F)) / 2;
A_Y_f = @(F) (F + jp1(F)) / 2;

lead1 = @(x) reshape(x, [1 size(x)]);
flat = @(F) F(:);

[nx, ny] = size(z);

%% Grid information
weight_u = OPT.RAWvec + zeros(OPT.nx, OPT.ny);
weight_v = OPT.RASvec + zeros(OPT.nx, OPT.ny);
dx = OPT.DXCvec;
dy = OPT.DYCsc;

%% Interpolation
% interpfn = @ppc_linterp;
interpfn = @ppc_pchip;
SppZ = interpfn(Z, S);
TppZ = interpfn(Z, T);
UppZ = interpfn(Z, U); % knowing Z is 1D vector, horizontally uniform
VppZ = interpfn(Z, V); % knowing Z is 1D vector, horizontally uniform

DATACUBE = OPT.data_cube;

%% slope error
tolz = OPT.tolz;
[sx, sy] = ntp_slope_error(SppZ, TppZ, Z, z, tolz, dx, dy, wrap);

% Calculate s2 and s
sx2 = sx .* sx .* (OPT.DXCvec .* OPT.DYGsc);
sy2 = sy .* sy .* (OPT.DXGvec .* OPT.DYCsc);
s2 = (sx2 + ip1(sx2) + sy2 + jp1(sy2)) ./ (2 * OPT.RACvec);
Df = 1000*s2;
s = sqrt(s2);

sx_rms = rms_Cgrid(0, sx, 0, 0, OPT.RAWvec, 0);
sy_rms = rms_Cgrid(0, 0, sy, 0, 0, OPT.RASvec);
s_rms  = rms_Cgrid(0, sx, sy, 0, OPT.RAWvec, OPT.RASvec);
s2_rms = rms_Cgrid(0, sx.^2, sy.^2, 0, OPT.RAWvec, OPT.RASvec);
Df_rms = rms_Cgrid(0, 1000*sx.^2, 1000*sy.^2, 0, OPT.RAWvec, OPT.RASvec);

if ~DATACUBE
    sx_st = sx(:,136:146);
    sy_st = sy(:,136:146);
    AreaX_st = OPT.RAWvec(:,136:146);
    AreaY_st = OPT.RASvec(:,136:146);
    sx_st_rms = rms_Cgrid(0, sx_st, 0, 0, AreaX_st, 0);
    sy_st_rms = rms_Cgrid(0, 0, sy_st, 0, 0, AreaY_st);
    s_st_rms  = rms_Cgrid(0, sx_st, sy_st, 0, AreaX_st, AreaY_st);
    s2_st_rms = rms_Cgrid(0, sx_st.^2, sy_st.^2, 0, AreaX_st, AreaY_st);
    Df_st_rms = rms_Cgrid(0, 1000*sx_st.^2, 1000*sy_st.^2, 0, AreaX_st, AreaY_st);
end

%% ehel
% u, v on the u and v grid
u = ppc_val(Z, UppZ, (z + im1(z)) / 2, 0);
v = ppc_val(Z, VppZ, (z + jm1(z)) / 2, 0);

u2 = u .* u .* (OPT.DXCvec .* OPT.DYGsc);
v2 = v .* v .* (OPT.DXGvec .* OPT.DYCsc);
u_t2 = (u2 + ip1(u2) + v2 + jp1(v2)) ./ (2 * OPT.RACvec);
u_t = sqrt(u_t2); % velocity on the tracer cell

% Calculate e^hel on (i-1/2, j):
e_hel_u = u.*sx;
% Calculate e^hel on (i, j-1/2):
e_hel_v = v.*sy;

e_hel_u_rms = rms_Cgrid(0, e_hel_u, 0, 0, OPT.RAWvec, 0);
e_hel_v_rms = rms_Cgrid(0, 0, e_hel_v, 0, 0, OPT.RASvec);

ehelx = u .* sx .* (OPT.DXCvec .* OPT.DYGsc);
ehely = v .* sy .* (OPT.DXGvec .* OPT.DYCsc);
e_hel = (ehelx + ip1(ehelx) + ehely + jp1(ehely)) ./ (2 * OPT.RACvec);

e_hel_rms = rms_Cgrid(e_hel, 0, 0, OPT.RACvec, 0, 0);

if ~DATACUBE
    e_hel_st = e_hel(:,136:146);
    Area_st = OPT.RACvec(:,136:146);
    
    e_hel_st_rms = rms_Cgrid(e_hel_st, 0, 0, Area_st, 0, 0);
end

%%
ad_df = 1e-3*Df + abs(e_hel);
ad_df_rms = rms_Cgrid_v2(e_hel, sx.^2, sy.^2, OPT.RACvec, OPT.RAWvec, OPT.RASvec);%method 1
% ad_df_rms = rms_Cgrid(ad_df, 0, 0, OPT.RACvec, 0, 0); % method 2
% ad_df_rms = rms_Cgrid(e_hel, sx.^2+ abs(e_hel), sy.^2+ abs(e_hel), OPT.RACvec, OPT.RAWvec, OPT.RASvec); %method 3

Tz = ppc_val_mex(Z, TppZ, z, 1);
ehelTz = Tz .* e_hel;
ehelTz_rms = rms_Cgrid(ehelTz, 0, 0, OPT.RACvec, 0, 0);

Sz = ppc_val_mex(Z, SppZ, z, 1);
ehelSz = Sz .* e_hel;
ehelSz_rms = rms_Cgrid(ehelSz, 0, 0, OPT.RACvec, 0, 0);

%% detail output for diagnostics
d.u = u;
d.v = v;
d.u_meanT = A_X_f(u);
d.v_meanT = A_Y_f(v);
d.u_t = u_t;
d.sx = sx;
d.sy = sy;
d.sx_rms = sx_rms;
d.sy_rms = sy_rms;
if ~DATACUBE
    d.s2_st_rms = s2_st_rms;
    d.e_hel_st_rms = e_hel_st_rms;
    d.Df_st_rms = Df_st_rms;
end
d.e_hel_u = e_hel_u;
d.e_hel_v = e_hel_v;
d.e_hel_u_rms = e_hel_u_rms;
d.e_hel_v_rms = e_hel_v_rms;
d.ehelTz = ehelTz;
d.ehelTz_rms = ehelTz_rms;
d.ehelSz = ehelSz;
d.ehelSz_rms = ehelSz_rms;
d.ad_df = ad_df;
d.ad_df_rms = ad_df_rms;
d.s2_rms = s2_rms;

end



