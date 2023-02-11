%% Load OCCA

% Edit `datafolder` as needed, to point to folder containing OCCA's *0406annclim..nc files
datafolder='/Users/yandonglang/OneDrive - UNSW/1. YL_PHD/3. Code/First year code/occa_data_calculation/OCCA/';
% datafolder='C:\Users\z5233717\OneDrive - UNSW\1. YL_PHD\3. Code\First year code\occa_data_calculation\OCCA\';

[g, SB, TB, pB, ETAN, ATMP, SAP] = load_OCCA(datafolder);

% Reshape the data into correct shape
tB = eos80_legacy_pt(SB,TB,zeros(size(SB)),pB); % 1980 code
SB = permute(SB, [3 1 2]);
tB = permute(tB, [3 1 2]); % in situ temperature
pB = permute(pB, [3 1 2]);
TB = permute(TB, [3 1 2]);
ZB = -g.RC; % ZB is depth and ZB >0

% get the arrray size
[nk, ni, nj] = size(SB);

% get the horizontal velocity
U = ncread(sprintf('%sDD%s.0406annclim.nc', datafolder, 'uvel'),'u');
U = permute(U, [3 1 2]); % depth x lon x lat
V = ncread(sprintf('%sDD%s.0406annclim.nc', datafolder, 'vvel'),'v');
V = permute(V, [3 1 2]); % depth x lon x lat

% interpolate salinity and temperature on the neutral density surface
% interpfn = @ppc_linterp;
interpfn = @ppc_pchip;
SppZ = interpfn(ZB, SB);
TppZ = interpfn(ZB, TB);

% SppP = interpfn(pB, SB);
% TppP = interpfn(pB, TB);


lead1 = @(x) reshape(x, [1 size(x)]);

im1 = @(F) circshift(F, [+1 0]);
jm1 = @(F) circshift(F, [0 +1]);
ip1 = @(F) circshift(F, [-1 0]);
jp1 = @(F) circshift(F, [0 -1]);

% avergece and backward difference
A_X = @(F) (F + im1(F)) / 2;
D_X_f = @(F)  ip1(F) - F;
A_Y = @(F) (F + jm1(F)) / 2;
D_Y_f = @(F)  jp1(F) - F;

A_X_f = @(F) (F + ip1(F)) / 2;
A_Y_f = @(F) (F + jp1(F)) / 2;

D_X_b = @(F)  F - im1(F);
D_Y_b = @(F)  F - jm1(F);

%% reference cast
% x0 = 243; y0 = -25; % South Pacific
x0 = 180; y0 = 0; % Pacific

[~,i0] = min(abs(g.XCvec - x0));
[~,j0] = min(abs(g.YCvec - y0));
I0 = sub2ind([g.nx,g.ny], i0, j0);

%% sigma surface
% Start with a potential density surface, crudely calculated
z_ref = 1200;
z0 = 2400;

SIGMA = eos(SB, TB, z_ref); % the potential density of the world 3D data reference to z_ref
sigma0 = interpfn(ZB, SIGMA(:,i0,j0), z0); % The isovalue of the potential density we want
z_sigma = interpfn(SIGMA, ZB, sigma0); % the potential density surface's depth of the potential density surface reference to z_ref with iso value sigma0
[s_sigma0, t_sigma0] = ppc_val2(ZB, SppZ, TppZ, z_sigma);

%% Neutral density labels

XC3 = repmat(permute(g.XCvec, [3 1 2]),[g.nz,1,g.ny]);
YC3 = repmat(permute(g.YCvec, [3 1 2]),[g.nz,g.nx,1]);
longOCCAzxy = XC3;
latOCCAzxy = YC3;
dpdz_bsq = (g.grav * g.rho_c * 1e-4);
pBB = dpdz_bsq * ZB;  % Get Boussinesq pressure (third argument to EOS) from Depth.

tic
% gn = gamma_n_2020(SB, tB, pB, longOCCAzxy, latOCCAzxy);
gscvE = gamma_scv(SB, tB, pBB, longOCCAzxy, latOCCAzxy); % ~ 3 mins
toc

%% Neutral density surface
g_iso = interpfn(ZB, gscvE(:,i0,j0), z0); % The isovalue of the potential density we want
% g_iso = 27.76;
[sscv, tscv, pscv] = eos80_legacy_neutral_surfaces(SB, tB, pBB, gscvE, g_iso);
tscv = eos80_legacy_pt(sscv,tscv,pscv, zeros(size(pscv))); % 1980 code, transfer in-situ temperature to potential temperature
sscv = squeeze(sscv);tscv = squeeze(tscv);pscv = squeeze(pscv);
zscv = pscv / dpdz_bsq; % transfer pressure to depth

%% gammaT
if exist('gsw_rho_CT_exact', 'file')
    SR = gsw_SR_from_SP(SB);
    SA = gsw_SA_from_SP(SB,pB,longOCCAzxy,latOCCAzxy);
    CT = gsw_CT_from_pt(SA,TB);
    tic
    [z_gt,p_gt,sigref,gammat] = gsw_gammat_analytic_os_2021(SR,CT); % 30 sec
    toc

    % gammaT surface
    gammat0 = interpfn(ZB, gammat(:,i0,j0), z0); % The isovalue of the potential density we want
    z_gammat = interpfn(gammat, ZB, gammat0); % the depth of the gammaT=27.5 surface, ZB is the depth
    [s_gammat, t_gammat] = ppc_val2(ZB, SppZ, TppZ, z_gammat);
else
    s_gammat = nan(ni, nj);
    t_gammat = nan(ni, nj);
    z_gammat = nan(ni, nj);
end

%% topobaric surface
OPTS = [];
[z_topo, s_topo, t_topo, RG, s0, t0, d_fn, diags_topo] = topobaric_surface(SB, TB, ZB, z_sigma, I0,[1;1], OPTS);

%% orthobaric surface
OPTS = [];
OPTS.REEB = false;

[z_ortho, s_ortho, t_ortho, RGortho, s0ortho, t0ortho, d_fnortho, diags_ortho] = topobaric_surface(SB, TB, ZB, z_sigma, I0,[1;1], OPTS);

%% omega_1
OPTS = [];
OPTS.POISSON = 0;
OPTS.INTERPFN = interpfn;

[zns, sns, tns, dns] = omega_surface(SB, TB, ZB, z_sigma, I0,[1;1], OPTS); %in x-y direction


%% Codegen for ppc_val* and ntp_slope
[nz, nx, ny] = size(SB);
t_Z = coder.typeof(0, [nz, 1], [0, 0]);
t_Sppc = coder.typeof(0, [4, nz-1, nx, ny], [1, 0, 0, 0]);
t_s = coder.typeof(0, [nx, ny], [1, 1]);
codegen('ppc_val', '-report', '-args', {t_Z, t_Sppc, t_s, 0})
codegen('ppc_val2', '-report', '-args', {t_Z, t_Sppc, t_Sppc, t_s, 0})
ni_ = max(4096, ni);
nj_ = max(4096, nj);
nk_ = max(128, nk);
ntp_slope_codegen(nk_, ni_, nj_, isvector(SB));

%% Grid information
OPTS.nx = g.nx; 
OPTS.ny = g.ny;
OPTS.DXCvec = g.DXCvec;
OPTS.DXGvec = g.DXGvec;
OPTS.DYCsc = g.DYCsc;
OPTS.DYGsc = g.DYGsc;
OPTS.RACvec = g.RACvec;
%OPTS.RACvec = 1;
OPTS.XCvec = g.XCvec;
OPTS.YCvec = g.YCvec;

% OPTS set up
OPTS.DAMP = 1;
% OPTS.x0 = 243; OPTS.y0 = -25; %reference cast in South Pacific
OPTS.x0 = 180; OPTS.y0 = 0; % Pacific
OPTS.FIGS_SHOW = 0;
OPTS.INTERPFN = @ppc_pchip;
OPTS.WRAP = [1;1];
OPTS.ITER = 50; % number of iterations

%% omega_hel
fprintf('omega_hel\n');
OPTS.data_cube=0;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;
OPTS.PCG = 0;
OPTS.LM = 0;
OPTS.H_SIM = 0;
OPTS.TK = 8e-9;

OPTS.CHOLESKY = 0;
OPTS.LSQLIN = 0;
OPTS.MODE = 1;

[z_hel,s_hel,t_hel, ~,d_hel] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% omega_hel_Tz
fprintf('omega_hel_Tz\n');
OPTS.data_cube=0;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;
OPTS.LM = 0;
OPTS.TK = 1e-11;

OPTS.CHOLESKY = 0;
OPTS.LSQLIN = 0;
OPTS.MODE = 9;

[z_helTz,s_helTz,t_helTz, ~,d_helTz] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% omega_hel_Sz
fprintf('omega_hel_Sz\n');
OPTS.data_cube=0;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;
OPTS.LM = 0;
OPTS.TK = 1e-13;

OPTS.CHOLESKY = 0;
OPTS.LSQLIN = 0;
OPTS.MODE = 10;

[z_helSz,s_helSz,t_helSz, ~,d_helSz] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);
           
%% omega_s
fprintf('omega_s\n');
OPTS.data_cube=0;
OPTS.STRAT = 1;
OPTS.SHEAR = 1;
% OPTS.LM = 1e-3;
OPTS.TK = 0;
OPTS.H_SIM = 0;
% OPTS.LM = 0;
% OPTS.TK = 1e-8;
OPTS.TOL_LSQR_REL = 1e-6;

OPTS.CHOLESKY = 0;
OPTS.MODE = 2;

[zns_s,sns_s,tns_s, ~, d_s] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% omega_s2xy
fprintf('omega_s2xy\n');
OPTS.STRAT = 1;
OPTS.SHEAR = 1;
OPTS.H_SIM = 0;
OPTS.data_cube=0;
OPTS.LM = 0;
% OPTS.TK = 1e-9;
OPTS.TK = 1e-9;

OPTS.CHOLESKY = 0;
OPTS.MODE = 6;

[zns_s2xy,sns_s2xy,tns_s2xy, ~, d_s2xy] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% omega_hel_s2xy 
fprintf('omega_hel_s2xy\n');
OPTS.data_cube=0;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;
OPTS.LM = 0;
OPTS.TK = 5e-9; % used for 20Aug2022 result
% OPTS.TK = 8e-9;
OPTS.H_SIM = 0;

% weight of ehelu and ehelv
OPTS.s2xyW = 1;
OPTS.FIGS_SHOW = 0;

OPTS.CHOLESKY = 0;
OPTS.LSQLIN = 0;
OPTS.MODE = 8;

[z_hels2xy,s_hels2xy,t_hels2xy, ~,d_hels2xy] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% save data
save omega_u_s_paper_06DEC_different_depth s_hel t_hel z_hel...
                               sns_s tns_s zns_s...
                               s_helTz t_helTz z_helTz...
                               s_helSz t_helSz z_helSz...
                               sns_s2xy tns_s2xy zns_s2xy...
                               s_hels2xy t_hels2xy z_hels2xy...
                               zns sns tns...
                               z_sigma s_sigma t_sigma...
                               tscv sscv zscv...
                               z_gammat s_gammat t_gammat
%     
