%% Data input
datafolder='/Users/yandonglang/OneDrive - UNSW/1. YL_PHD/3. Code/First year code/occa_data_calculation/OCCA/';
[g, SB, TB, ~, ETAN, ATMP, SAP] = load_OCCA(datafolder);

g_iso = 27.76;

% add something
% Reshape the data into correct shape
SB = permute(SB, [3 1 2]);
%pB = permute(pB, [3 1 2]);  % this is OCCA's dynamic pressure.  Not what actually gets used in the EOS
TB = permute(TB, [3 1 2]);
ZB = -g.RC; % ZB is depth and ZB >0

dpdz_bsq = (g.grav * g.rho_c * 1e-4);
pB = dpdz_bsq * ZB;  % Get Boussinesq pressure (third argument to EOS) from Depth.

tB = eos80_legacy_pt(SB,TB,zeros(size(SB)),pB);  % compute in-situ temperature from potential temperature

% SPPX, TppX
%interpfn = @ppc_pchip;
interpfn = @ppc_linterp;  % Better to use linear since this is what the gamma^n code assumes -- I think?
SppZ = interpfn(ZB, SB);
TppZ = interpfn(ZB, TB);

% get the arrray size
[nk, ni, nj] = size(SB);


if 0  % After you've done this once, can change to 0 so functions are not constantly re-compiled
  PATH_NS = '~/work/projects-gfd/neutral-surfaces/'; % << ADJUST AS NEEDED >>
  PATH_EOS = [PATH_NS '/lib/eos/active/'];
  eoscg_set_bsq_param([PATH_NS 'lib/eos/eoscg_densjmd95_bsq.m'   ] , [PATH_EOS 'eos.m'  ], g.grav, g.rho_c);
  eoscg_set_bsq_param([PATH_NS 'lib/eos/eoscg_densjmd95_bsq_dz.m'] , [PATH_EOS 'eos_p.m'], g.grav, g.rho_c);  % not needed for this code, but here it is anyway
  eoscg_set_bsq_param([PATH_NS 'lib/eos/eoscg_densjmd95_bsq_s_t.m'], [PATH_EOS 'eos_s_t.m'], g.grav, g.rho_c);  % not needed for this code, but here it is anyway
end

assert(abs(eos(35,2,1000) - 1.032659160630248e+03) < 1e-12, 'Wrong eos.m in use!  Expected the Boussinesq version of densjmd95 with grav = %.8f and rho_c = %.8f', g.grav, g.rho_c)

%% Neutral density labels for each bottle

XC3 = repmat(permute(g.XCvec, [3 1 2]),[g.nz,1,g.ny]);
YC3 = repmat(permute(g.YCvec, [3 1 2]),[g.nz,g.nx,1]);
longOCCAzxy = XC3;
latOCCAzxy = YC3;
% gn = eos80_legacy_gamma_n(SB, tB, pB, longOCCAzxy, latOCCAzxy);  % ~ 2 mins
gn = eos80_legacy_gamma_n_STP_GJS_no_extrapolation(SB, tB, pB, longOCCAzxy, latOCCAzxy);
gscv = eos80_legacy_gamma_scv_STP_GJS_no_extrapolation(SB, tB, pB, longOCCAzxy, latOCCAzxy);

%% gamma^n surface
[s, t, p] = eos80_legacy_neutral_surfaces(SB, tB, pB, gscv, g_iso);
t = eos80_legacy_pt(s,t,p, zeros(size(p))); % 1980 code, transfer in-situ temperature to potential temperature
s = squeeze(s);
t = squeeze(t);
p = squeeze(p);
zscv = p / dpdz_bsq;  % convert Boussinesq pressure back to Depth

%% slope error VENM for gamma^n surface
[sx, sy] = ntp_slope_error(SppZ, TppZ, ZB, z, 1e-8, g.DXCvec, g.DYCsc, [1;0]);

nanrms(sx(:));  % 1.516352307954056e-05

%% gammaT
SR = gsw_SR_from_SP(SB);
pBzxy = repmat(pB, [1, g.nx, g.ny]);
SA = gsw_SA_from_SP(SB,pBzxy,longOCCAzxy,latOCCAzxy);
CT = gsw_CT_from_pt(SA,TB);
[z_gt,p_gt,sigref,gammat] = gsw_gammat_analytic_os_2021(SR,CT); % < 1 minute

%% gammaT surface
z_gammat = interpfn(gammat, ZB, g_iso); % the depth of the gammaT=27.5 surface, ZB is the depth
s_gammat = interpfn(ZB, SB, z_gammat);
t_gammat = interpfn(ZB, TB, z_gammat);

%% slope error VENM
z = z_gammat;
[sx, sy] = ntp_slope_error(SppZ, TppZ, ZB, z, 1e-8, g.DXCvec, g.DYCsc, [1;0]);

nanrms(sx(:)) % 1.518520612139688e-05
