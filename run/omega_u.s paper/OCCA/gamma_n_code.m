%% Data input

[g, SB, TB, pB, ETAN, ATMP, SAP] = load_OCCA(datafolder); % please specify the data folder

% add something
% Reshape the data into correct shape
tB = eos80_legacy_pt(SB,TB,zeros(size(SB)),pB); % 1980 code
SB = permute(SB, [3 1 2]);
tB = permute(tB, [3 1 2]); % in situ temperature
pB = permute(pB, [3 1 2]);
TB = permute(TB, [3 1 2]);
ZB = -g.RC; % ZB is depth and ZB >0

% SPPX, TppX
interpfn = @ppc_pchip;
SppZ = interpfn(ZB, SB);
TppZ = interpfn(ZB, TB);

% get the arrray size
[nk, ni, nj] = size(SB);

%% Neutral density labels

XC3 = repmat(permute(g.XCvec, [3 1 2]),[g.nz,1,g.ny]);
YC3 = repmat(permute(g.YCvec, [3 1 2]),[g.nz,g.nx,1]);
longOCCAzxy = XC3;
latOCCAzxy = YC3;
% label with gamma^n
[gscvE, ~, ~, SE, tE, TE, pE] = eos80_legacy_gamma_n_STP_GJS_no_extrapolation(SB, tB, pBB, longOCCAzxy, latOCCAzxy);

%% Neutral density surface 
g_iso = 27.76;
[sscv, tscv, pscv] = eos80_legacy_neutral_surfaces(SE, tE, pE, gscvE, g_iso);
tscv = eos80_legacy_pt(sscv,tscv,pscv, zeros(size(pscv))); % 1980 code, transfer in-situ temperature to potential temperature
sscv = squeeze(sscv);tscv = squeeze(tscv);pscv = squeeze(pscv);
zscv = -gsw_z_from_p(pscv,g.YCvec); % transfer pressure to depth

%% slope error VENM
[sx, sy] = ntp_slope_error(SppZ, TppZ, ZB, zscv, 1e-8, g.DXCvec, g.DYCsc, [1;0]);

nanrms(sx(:))

%% gammaT
SR = gsw_SR_from_SP(SB);
SA = gsw_SA_from_SP(SB,pB,longOCCAzxy,latOCCAzxy);
CT = gsw_CT_from_pt(SA,TB);
[z_gt,p_gt,sigref,gammat] = gsw_gammat_analytic_os_2021(SR,CT);

%% gammaT surface
z_gammat = interpfn(gammat, ZB, 27.76); % the depth of the gammaT=27.5 surface, ZB is the depth
s_gammat = interpfn(ZB, SB, z_gammat);
t_gammat = interpfn(ZB, TB, z_gammat);

%%
% for gamma^T
[ex,ey,sx,sy] = ntp_errors(s_gammat, t_gammat, z_gammat, g.DXCvec, g.DYCsc,1,0,[1;0], {},9.81, SB, TB,ZB);
% for gamma^n
[ex1,ey1,sx1,sy1] = ntp_errors(sscv, tscv, zscv, g.DXCvec, g.DYCsc,1,0,[1;0], {},9.81, SB, TB,ZB);

figure('Position', [0 0 1000 1000])
subplot(2,1,1)
pcolor(ex')
shading flat
colorbar
caxis([-1e-7, 1e-7])
title('$e_x$ of $\gamma^T$', 'Interpreter', 'latex')
set(gca,'fontsize', 16);

subplot(2,1,2)
pcolor(ex1')
shading flat
colorbar
caxis([-1e-7, 1e-7])
title('$e_x$ of $\gamma^n$', 'Interpreter', 'latex')
set(gca,'fontsize', 16);
