%% Load OCCA
datafolder='/Users/yandonglang/OneDrive - UNSW/1. YL_PHD/3. Code/First year code/occa_data_calculation/OCCA/';
% datafolder='C:\Users\z5233717\OneDrive - UNSW\1. YL_PHD\3. Code\First year code\occa_data_calculation\OCCA\';

[g, SB, TB, pB, ETAN, ATMP, SAP] = load_OCCA(datafolder);

% add something
% Reshape the data into correct shape
tB = eos80_legacy_pt(SB,TB,zeros(size(SB)),pB); % 1980 code
SB = permute(SB, [3 1 2]);
tB = permute(tB, [3 1 2]);
pB = permute(pB, [3 1 2]);
TB = permute(TB, [3 1 2]);
ZB = -g.RC; % ZB is depth and ZB >0

% get the arrray size
[nk, ni, nj] = size(SB);

% get the horizontal velocity
U = ncread('DDuvel.0406clim.nc','u');U = squeeze(U(:,:,:,1));U = permute(U, [3 1 2]);
V = ncread('DDvvel.0406clim.nc','v');V = squeeze(V(:,:,:,1));V = permute(V, [3 1 2]);

% interpolate salinity and temperature on the neutral density surface
% interpfn = @ppc_linterp;
interpfn = @ppc_pchip;
SppZ = interpfn(ZB, SB);
TppZ = interpfn(ZB, TB);

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
z_ref = 750;
z0 = 1500;

SIGMA = eos(SB, TB, z_ref); % the potential density of the world 3D data reference to z_ref
sigma0 = interpfn(ZB, SIGMA(:,i0,j0), z0); % The isovalue of the potential density we want
z_sigma = interpfn(SIGMA, ZB, sigma0); % the potential density surface's depth of the potential density surface reference to z_ref with iso value sigma0
s_sigma = interpfn(ZB, SB, z_sigma);
t_sigma = interpfn(ZB, TB, z_sigma);

z_ref = 500;

SIGMA = eos(SB, TB, z_ref); % the potential density of the world 3D data reference to z_ref
sigma0 = interpfn(ZB, SIGMA(:,i0,j0), z0); % The isovalue of the potential density we want
z_sigma500 = interpfn(SIGMA, ZB, sigma0); % the potential density surface's depth of the potential density surface reference to z_ref with iso value sigma0
s_sigma500 = interpfn(ZB, SB, z_sigma500);
t_sigma500 = interpfn(ZB, TB, z_sigma500);

z_ref = 1000;

SIGMA = eos(SB, TB, z_ref); % the potential density of the world 3D data reference to z_ref
sigma0 = interpfn(ZB, SIGMA(:,i0,j0), z0); % The isovalue of the potential density we want
z_sigma1000 = interpfn(SIGMA, ZB, sigma0); % the potential density surface's depth of the potential density surface reference to z_ref with iso value sigma0
s_sigma1000 = interpfn(ZB, SB, z_sigma1000);
t_sigma1000 = interpfn(ZB, TB, z_sigma1000);

z_ref = 2000;

SIGMA = eos(SB, TB, z_ref); % the potential density of the world 3D data reference to z_ref
sigma0 = interpfn(ZB, SIGMA(:,i0,j0), z0); % The isovalue of the potential density we want
z_sigma2000 = interpfn(SIGMA, ZB, sigma0); % the potential density surface's depth of the potential density surface reference to z_ref with iso value sigma0
s_sigma2000 = interpfn(ZB, SB, z_sigma2000);
t_sigma2000 = interpfn(ZB, TB, z_sigma2000);

%% In-situ density anomaly
z0 = 2000; % the depth that I want the iso-value of the in-situ density surface
insitu_density = eos(SB, TB, ZB); % caculate the ocean's in situ density
insitu0 = interpfn(ZB, insitu_density(:,i0,j0), z0); % The isovalue of the in situ density surface we want, i0, j0 are the reference cast
z_insitu = interpfn(insitu_density, ZB, insitu0); % the in situ density surface depth

%% topobaric surface
OPTS = [];
[z_topo, s_topo, t_topo, RG, s0, t0, d_fn, diags_topo] = topobaric_surface(SB, TB, ZB, z_sigma, I0,[1;1], OPTS);

%% omega_1
OPTS = [];

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
OPTS.XCvec = g.XCvec;
OPTS.YCvec = g.YCvec;

% OPTS set up
OPTS.DAMP = 1;
% OPTS.x0 = 243; OPTS.y0 = -25; %reference cast in South Pacific
OPTS.x0 = 180; OPTS.y0 = 0; % Pacific
OPTS.FIGS_SHOW = 0;
OPTS.INTERPFN = @ppc_pchip;
OPTS.WRAP = [1;1];

%% omega_hel
OPTS.data_cube=0;
OPTS.ITER = 100;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;

OPTS.LM = 0;
OPTS.TK = 1e-5;
%OPTS.TK = 0;

% OPTS.TK = 0;
% OPTS.LM = 1;

OPTS.MODE = 1;

[z_hel,s_hel,t_hel, ~,d_hel] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% omega_hel_Tz
OPTS.data_cube=0;
OPTS.ITER = 100;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;

OPTS.LM = 0;
OPTS.TK = 1e-9;

OPTS.MODE = 9;

[z_helTz,s_helTz,t_helTz, ~,d_helTz] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% omega_hel_Sz
OPTS.data_cube=0;
OPTS.ITER = 100;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;

OPTS.LM = 0;
OPTS.TK = 1e-10;

OPTS.MODE = 10;

[z_helSz,s_helSz,t_helSz, ~,d_helSz] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);
           
%% omega_s
OPTS.MODE = 2;
OPTS.STRAT = 1;
OPTS.SHEAR = 1;
% OPTS.LM = 1e-3;
OPTS.TK = 0;
OPTS.H_SIM = 1;
% OPTS.LM = 0;
% OPTS.TK = 1e-8;

OPTS.ITER = 100;

[zns_s,sns_s,tns_s, ~, d_s] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% omega_s2xy
OPTS.MODE = 6;
OPTS.STRAT = 1;
OPTS.SHEAR = 1;
OPTS.ITER = 100;
OPTS.H_SIM = 0;

OPTS.LM = 0;
OPTS.TK = 1e-9;

% OPTS.TK = 0;
% OPTS.LM = 1;

[zns_s2xy,sns_s2xy,tns_s2xy, ~, d_s2xy] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% omega_hel_s2xy 
OPTS.ITER = 100;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;

OPTS.LM = 0;
OPTS.TK = 1e-5;
OPTS.H_SIM = 0;

% OPTS.TK = 0;
% OPTS.LM = 1;

% weight of ehelu and ehelv
OPTS.s2xyW = 1;
OPTS.FIGS_SHOW = 0;

OPTS.MODE = 8;

[z_hels2xy,s_hels2xy,t_hels2xy, ~,d_hels2xy] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% s, e^hel and flux  
% prepare for the grid information
OPT.DXCvec = g.DXCvec; % dx
OPT.DYCsc = g.DYCsc; % dy
OPT.DYGsc = g.DYGsc;
OPT.DXGvec = g.DXGvec;
OPT.RAWvec = g.RAWvec; % area of the grid in u grid
OPT.RASvec = g.RASvec; % area of the grid in v grid
OPT.RACvec = g.RACvec; % area of the grid in tracer grid
OPT.nx = g.nx;
OPT.ny = g.ny;
OPT.tolz = 1e-8; % tolerance of slope error

% s, ehel and flux for omega+
OPT.MODES = 0; % using epsilon to calculate slope errors
OPT.MODEHEL = 1; % Using divergence method to calculate e^hel
[s_ns, s2ns, sns_rms, s2ns_rms, ehelns, ehelns_rms, fluxns, detailns] =...
                                 s_ehel_flux(sns,tns, zns, SB, TB, ZB, U, V, [1;1], OPT);

% s, ehel and flux for omega_hel
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % using divergence method to calculate e^hel
[ss_hel, s2hel, ss_hel_rms, s2hel_rms, ehel_hel, ehel_hel_rms, flux_hel, detail_hel] =...
                            s_ehel_flux(s_hel,t_hel, z_hel, SB, TB, ZB, U, V, [1;1], OPT);

% s, ehel and flux for omega_helTz
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % using divergence method to calculate e^hel
[ss_helTz, s2helTz, ss_helTz_rms, s2helTz_rms, ehel_helTz, ehel_helTz_rms, flux_helTz, detail_helTz] =...
                            s_ehel_flux(s_helTz,t_helTz, z_helTz, SB, TB, ZB, U, V, [1;1], OPT);

% s, ehel and flux for omega_helSz
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % using divergence method to calculate e^hel
[ss_helSz, s2helSz, ss_helSz_rms, s2helSz_rms, ehel_helSz, ehel_helSz_rms, flux_helSz, detail_helSz] =...
                            s_ehel_flux(s_helSz,t_helSz, z_helSz, SB, TB, ZB, U, V, [1;1], OPT);                        
                                             
% s, ehel and flux for omega_s
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % using divergence method to calculate e^hel
[s_ns_s, s2ns_s, sns_s_rms, s2ns_s_rms, ehelns_s, ehelns_s_rms, fluxns_s, detailns_s] =...
                             s_ehel_flux(sns_s,tns_s, zns_s, SB, TB, ZB, U, V, [1;1], OPT);
                         
% s, ehel and flux for topobaric surface
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % using divergence method to calculate e^hel
[stopo, s2_topo, stopo_rms, s2_topo_rms, ehel_topo, ehel_topo_rms, flux_topo, detail_topo] =...
                             s_ehel_flux(s_topo,t_topo, z_topo, SB, TB, ZB, U, V, [1;1], OPT);
                         

% s, ehel and flux for sigma
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % Using divergence method to calculate e^hel
[s_pt, s2pt, spt_rms, s2pt_rms, ehelpt, ehelpt_rms, fluxpt, detailpt] =...
                                 s_ehel_flux(s_sigma,t_sigma, z_sigma, SB, TB, ZB, U, V, [1;1], OPT);
                             
% s, ehel and flux for sigma500
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % Using divergence method to calculate e^hel
[s_pt5, s2pt5, spt5_rms, s2pt5_rms, ehelpt5, ehelpt5_rms, fluxpt5, detailpt5] =...
                                 s_ehel_flux(s_sigma500,t_sigma500, z_sigma500, SB, TB, ZB, U, V, [1;1], OPT);                             

% s, ehel and flux for sigma1000
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % Using divergence method to calculate e^hel
[s_pt10, s2pt10, spt10_rms, s2pt10_rms, ehelpt10, ehelpt10_rms, fluxpt10, detailpt10] =...
                                 s_ehel_flux(s_sigma1000,t_sigma1000, z_sigma1000, SB, TB, ZB, U, V, [1;1], OPT);                             

% s, ehel and flux for sigma1000
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % Using divergence method to calculate e^hel
[s_pt20, s2pt20, spt20_rms, s2pt20_rms, ehelpt20, ehelpt20_rms, fluxpt20, detailpt20] =...
                                 s_ehel_flux(s_sigma2000, t_sigma2000, z_sigma2000, SB, TB, ZB, U, V, [1;1], OPT);                             
                              
% s, ehel and flux for omega_s^2 seperate in u and v
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % using divergence method to calculate e^hel
[s_ns_s2xy, s2ns_s2xy, sns_s2xy_rms, s2ns_s2xy_rms, ehelns_s2xy, ehelns_s2xy_rms, fluxns_s2xy, detailns_s2xy] =...
                             s_ehel_flux(sns_s2xy,tns_s2xy, zns_s2xy, SB, TB, ZB, U, V, [1;1], OPT);                          
  
% s, ehel and flux for omega_hel_s2uv
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % using divergence method to calculate e^hel
[ss_hel_s2xy, s2hel_s2xy, ss_hel_s2xy_rms, s2hel_s2xy_rms, ehel_hel_s2xy, ehel_hel_s2xy_rms, flux_hel_s2xy, detail_hel_s2xy] =...
                            s_ehel_flux(s_hels2xy,t_hels2xy, z_hels2xy, SB, TB, ZB, U, V, [1;1], OPT);  
                        
% make a table
surfaces  = {'omega+', 'omega_hel', 'omega_s', 'omega_s2xy',  'omega_hel_s2xy'}';
slope_rms = [sns_rms;   ss_hel_rms;  sns_s_rms;   sns_s2xy_rms;   ss_hel_s2xy_rms];
s2_rms    = [s2ns_rms;  s2hel_rms;   s2ns_s_rms;    s2ns_s2xy_rms;    s2hel_s2xy_rms];
fluxes    = [fluxns;    flux_hel;    fluxns_s;     fluxns_s2xy;     flux_hel_s2xy];
ehel_rms  = [ehelns_rms;ehel_hel_rms;ehelns_s_rms; ehelns_s2xy_rms;ehel_hel_s2xy_rms];

T = table(surfaces,slope_rms,s2_rms,ehel_rms, fluxes)  

surfacess  = {'sigma', 'omega+', 'omega_s', 'omega_hel', 'omega_helTz', 'omega_helSz', 'omega_s2xy',  'omega_hel_s2xy'}';
ehel_south_rms = [detailpt.e_hel_st_rms; detailns.e_hel_st_rms; detailns_s.e_hel_st_rms; detail_hel.e_hel_st_rms; detail_helTz.e_hel_st_rms; detail_helSz.e_hel_st_rms; detailns_s2xy.e_hel_st_rms; detail_hel_s2xy.e_hel_st_rms];
Df_south_rms =   [detailpt.s2_st_rms; detailns.s2_st_rms; detailns_s.s2_st_rms; detail_hel.s2_st_rms; detail_helTz.s2_st_rms; detail_helSz.s2_st_rms; detailns_s2xy.s2_st_rms; detail_hel_s2xy.s2_st_rms];

TT = table(surfacess,ehel_south_rms,Df_south_rms)  

%% Figure 1:e^hel maps
OPTS_FIGS.LANDCOL = [1 1 1]*0; % Black
OPTS_FIGS.NANCOL = [1 1 1]* .75; % Grey
OPTS_FIGS.LATLIM = [-80, 72]; % The z(180,0) = 2000 surfaces exclude the Arctic.
land_OCCA = squeeze(isnan(SB(1,:,:)));

OPTS_plot.nx = 4;
OPTS_plot.ny = 2;
OPTS_plot.X = g.XCvec;
OPTS_plot.Y = g.YCvec;
OPTS_plot.land_mask = land_OCCA;
OPTS_plot.bwr = ones(1,8);
OPTS_plot.font_size = 15;
OPTS_plot.AXES = {'Margin', .04, 'Spacing', .03};

figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -1e-7*ones(1,8);
OPTS_plot.p_lim = 1e-7*ones(1,8);

data = cat(3, ehelpt, ehelns, ehelns_s, ehel_hel, ehel_helTz, ehel_helSz, ehelns_s2xy, ehel_hel_s2xy);
data_rms = [ehelpt_rms, ehelns_rms, ehelns_s_rms, ehel_hel_rms, ehel_helTz_rms,...
                                                       ehel_helSz_rms, ehelns_s2xy_rms, ehel_hel_s2xy_rms];

quantity_txt = '$e^{hel}$ on ';        
rms_txt = 'RMS of $e^{hel}$ = %.2d m/s';        
txt1 = ['(a) ', quantity_txt,'$\\sigma$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{hel}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{hel\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{helS_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{s^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{hel+s^2}$-surface, ', rms_txt];

title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 1:e^hel maps
OPTS_FIGS.LANDCOL = [1 1 1]*0; % Black
OPTS_FIGS.NANCOL = [1 1 1]* .75; % Grey
OPTS_FIGS.LATLIM = [-80, 72]; % The z(180,0) = 2000 surfaces exclude the Arctic.
land_OCCA = squeeze(isnan(SB(1,:,:)));

OPTS_plot.nx = 4;
OPTS_plot.ny = 2;
OPTS_plot.X = g.XCvec;
OPTS_plot.Y = g.YCvec;
OPTS_plot.land_mask = land_OCCA;
OPTS_plot.bwr = ones(1,8);
OPTS_plot.font_size = 15;
OPTS_plot.AXES = {'Margin', .04, 'Spacing', .03};

figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -1e-7*ones(1,8);
OPTS_plot.p_lim = 1e-7*ones(1,8);

data = cat(3, ehel_topo, ehelns, ehelns_s, ehel_hel, ehel_helTz, ehel_helSz, ehelns_s2xy, ehel_hel_s2xy);
data_rms = [ehel_topo_rms, ehelns_rms, ehelns_s_rms, ehel_hel_rms, ehel_helTz_rms,...
                                                       ehel_helSz_rms, ehelns_s2xy_rms, ehel_hel_s2xy_rms];

quantity_txt = '$e^{hel}$ on ';        
rms_txt = 'RMS of $e^{hel}$ = %.2d m/s';        
txt1 = ['(a) ', quantity_txt,'$\\tau$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{hel}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{hel\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{helS_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{s^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{hel+s^2}$-surface, ', rms_txt];

title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);
                            
%% e_hel Tz map
Tppc = ppc_pchip(ZB,TB);
Tz_ns = ppc_val_mex(ZB, Tppc, zns, 1);
Tz_pt = ppc_val_mex(ZB, Tppc, z_sigma, 1);
Tz_ns_s = ppc_val_mex(ZB, Tppc, zns_s, 1);
Tz_helTz = ppc_val_mex(ZB, Tppc, z_helTz, 1);
Tz_helSz = ppc_val_mex(ZB, Tppc, z_helSz, 1);
Tz_hel = ppc_val_mex(ZB, Tppc, z_hel, 1);
Tz_hels2xy = ppc_val_mex(ZB, Tppc, z_hels2xy, 1);
Tz_s2xy = ppc_val_mex(ZB, Tppc, zns_s2xy, 1);

ehelTz_pt = Tz_pt .* ehelpt;
ehelTz_ns = Tz_ns .* ehelns;
ehelTz_ns_s = Tz_ns_s .* ehelns_s;
ehelTz_helTz = Tz_helTz .* ehel_helTz;
ehelTz_helSz = Tz_helSz .* ehel_helSz;
ehelTz_hels2xy = Tz_hels2xy .* ehel_hel_s2xy;
ehelTz_s2xy = Tz_s2xy .* ehelns_s2xy;
ehelTz_hel = Tz_hel .* ehel_hel;

figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -1e-9*ones(1,8);
OPTS_plot.p_lim = 1e-9*ones(1,8);
OPTS_plot.bwr = ones(1,8);

data = cat(3, ehelTz_pt, ehelTz_ns, ehelTz_ns_s, ehelTz_hel,...
                                                 ehelTz_helTz, ehelTz_helSz, ehelTz_s2xy, ehelTz_hels2xy);
data_rms = [root_mean_square(ehelTz_pt(:)),    root_mean_square(ehelTz_ns(:)),...
            root_mean_square(ehelTz_ns_s(:)),  root_mean_square(ehelTz_hel(:)),...
            root_mean_square(ehelTz_helTz(:)), root_mean_square(ehelTz_helSz(:)),...
            root_mean_square(ehelTz_s2xy(:)),  root_mean_square(ehelTz_hels2xy(:))];
quantity_txt = '$\\Theta_z e^{hel}$ on ';        
rms_txt = 'RMS of $\\Theta_z e^{hel}$ = %.2d $^\\circ$C/s';        
txt1 = ['(a) ', quantity_txt,'$\\sigma$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{hel}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{hel\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{helS_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{s^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{hel+s^2}$-surface, ', rms_txt];
title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% e_hel Sz map
Sppc = ppc_pchip(ZB,SB);
Sz_ns = ppc_val_mex(ZB, Sppc, zns, 1);
Sz_pt = ppc_val_mex(ZB, Sppc, z_sigma, 1);
Sz_ns_s = ppc_val_mex(ZB, Sppc, zns_s, 1);
Sz_helTz = ppc_val_mex(ZB, Sppc, z_helTz, 1);
Sz_helSz = ppc_val_mex(ZB, Sppc, z_helSz, 1);
Sz_hel = ppc_val_mex(ZB, Sppc, z_hel, 1);
Sz_hels2xy = ppc_val_mex(ZB, Sppc, z_hels2xy, 1);
Sz_s2xy = ppc_val_mex(ZB, Sppc, zns_s2xy, 1);

ehelSz_pt = Sz_pt .* ehelpt;
ehelSz_ns = Sz_ns .* ehelns;
ehelSz_ns_s = Sz_ns_s .* ehelns_s;
ehelSz_helTz = Sz_helTz .* ehel_helTz;
ehelSz_helSz = Sz_helSz .* ehel_helSz;
ehelSz_hels2xy = Sz_hels2xy .* ehel_hel_s2xy;
ehelSz_s2xy = Sz_s2xy .* ehelns_s2xy;
ehelSz_hel = Sz_hel .* ehel_hel;

figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -1e-10*ones(1,8);
OPTS_plot.p_lim = 1e-10*ones(1,8);
OPTS_plot.bwr = ones(1,8);

data = cat(3, ehelSz_pt, ehelSz_ns, ehelSz_ns_s, ehelSz_hel,...
                                                 ehelSz_helTz, ehelSz_helSz, ehelSz_s2xy, ehelSz_hels2xy);
data_rms = [root_mean_square(ehelSz_pt(:)),    root_mean_square(ehelSz_ns(:)),...
            root_mean_square(ehelSz_ns_s(:)),  root_mean_square(ehelSz_hel(:)),...
            root_mean_square(ehelSz_helTz(:)), root_mean_square(ehelSz_helSz(:)),...
            root_mean_square(ehelSz_s2xy(:)),  root_mean_square(ehelSz_hels2xy(:))];
quantity_txt = '$S_z e^{hel}$ on ';        
rms_txt = 'RMS of $S_z e^{hel}$ = %.2d psu/s';        
txt1 = ['(a) ', quantity_txt,'$\\sigma$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{hel}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{hel\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{helS_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{s^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{hel+s^2}$-surface, ', rms_txt];
title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 2: Slope Error map
figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -8*ones(1,8);
OPTS_plot.p_lim = -3*ones(1,8);
OPTS_plot.bwr = zeros(1,8);

data = cat(3, log10(s_pt), log10(s_ns), log10(s_ns_s), log10(ss_hel),...
                                   log10(ss_helTz), log10(ss_helSz), log10(s_ns_s2xy), log10(ss_hel_s2xy));
data_rms = [spt_rms, sns_rms, sns_s_rms, ss_hel_rms, ss_helTz_rms,...
                                                       ss_helSz_rms, sns_s2xy_rms, ss_hel_s2xy_rms];
quantity_txt = '$log_{10} \\bf(s)$ on ';        
rms_txt = 'RMS = %.2d';        
txt1 = ['(a) ', quantity_txt,'$\\sigma$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{hel}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{hel\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{helS_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{s^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{hel+s^2}$-surface, ', rms_txt];
title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);


%% Figure 2: Slope Error map
figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -8*ones(1,8);
OPTS_plot.p_lim = -3*ones(1,8);
OPTS_plot.bwr = zeros(1,8);

data = cat(3, log10(stopo), log10(s_ns), log10(s_ns_s), log10(ss_hel),...
                                   log10(ss_helTz), log10(ss_helSz), log10(s_ns_s2xy), log10(ss_hel_s2xy));
data_rms = [stopo_rms, sns_rms, sns_s_rms, ss_hel_rms, ss_helTz_rms,...
                                                       ss_helSz_rms, sns_s2xy_rms, ss_hel_s2xy_rms];
quantity_txt = '$log_{10} \\bf(s)$ on ';        
rms_txt = 'RMS = %.2d';        
txt1 = ['(a) ', quantity_txt,'$\\tau$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{hel}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{hel\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{helS_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{s^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{hel+s^2}$-surface, ', rms_txt];
title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 3: K s^2 map
df_ns = 1000*s2ns;
df_ns_s = 1000*s2ns_s;
df_pt = 1000*s2pt;
df_hel = 1000*s2hel;
df_helTz = 1000*s2helTz;
df_helSz = 1000*s2helSz;
df_s2xy = 1000*s2ns_s2xy;
df_hel_s2xy = 1000*s2hel_s2xy;

figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -16*ones(1,8);
OPTS_plot.p_lim = -5*ones(1,8);
OPTS_plot.bwr = zeros(1,8);

data = cat(3, log10(df_pt), log10(df_ns), log10(df_ns_s), log10(df_hel),...
                                   log10(df_helTz), log10(df_helSz), log10(df_s2xy), log10(df_hel_s2xy));
data_rms = [1000*s2pt_rms, 1000*s2ns_rms, 1000*s2ns_s_rms, 1000*s2hel_rms, 1000*s2helTz_rms,...
                                              1000*s2helSz_rms, 1000*s2ns_s2xy_rms, 1000*s2hel_s2xy_rms];
quantity_txt = '$log_{10} (D^f)$ on ';        
rms_txt = 'RMS = %.2d m$^2$ s$^{-1}$';        
txt1 = ['(a) ', quantity_txt,'$\\sigma$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{hel}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{hel\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{helS_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{s^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{hel+s^2}$-surface, ', rms_txt];
title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 4: The plot of cos of omega+, omega_hel and sigma
sintheta_ns = ehelns./(abs(s_ns).*abs(detailns.u_t));
sintheta_hel = ehel_hel./(abs(ss_hel).*abs(detail_hel.u_t));
sintheta_pt = ehelpt./(abs(s_pt).*abs(detailpt.u_t));


OPTS_AXES = {'Margin', .06, 'Spacing', .03};
figure('Position', [0 0 800 1800])

ax1 = subaxis(3,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, sintheta_ns , land_OCCA, OPTS_FIGS);
colorbar(ax1)
txt = '(a) $v.s/|v||s|$  on $\\omega_+$-surface';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
caxis([-1, 1])
ax1.FontSize = 15;
ax1.XTickLabel = [];

ax2 = subaxis(3,1,2, OPTS_AXES{:});
hf = fig_map(ax2, g.XCvec, g.YCvec, sintheta_hel , land_OCCA, OPTS_FIGS);
colorbar(ax2)
txt = '(b) $v.s/|v||s|$  on $\\omega_{hel}$-surface';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
caxis([-1, 1])
cm = parula(64);
cm = [cm(1:end-1,:); flipud(cm)];
colormap(cm);
ax2.FontSize = 15;
ax2.XTickLabel = []; 

ax3 = subaxis(3,1,3, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, sintheta_pt , land_OCCA, OPTS_FIGS);
colorbar(ax3)
txt = '(c) $v.s/|v||s|$  on $\\sigma$-surface';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
caxis([-1, 1])
cm = parula(64);
cm = [cm(1:end-1,:); flipud(cm)];
colormap(cm);
ax3.FontSize = 15;

%% Figure 5:Depth map
figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = [0,   -50, -10, -10, -5, -5, -5, -5];
OPTS_plot.p_lim = [1800, 50,  10,  10,  5,  5,  5,  5];
OPTS_plot.bwr = [0, 1, 1, 1, 1, 1, 1, 1];

data = cat(3, zns, (z_sigma-zns), (zns_s-zns), (z_hel-zns),...
                                       (z_helTz-z_hel), (z_helSz-z_hel), (zns_s2xy-zns), (z_hels2xy-z_hel));
data_rms = [root_mean_square(zns),          root_mean_square(z_sigma-zns),...
            root_mean_square(zns_s-zns),    root_mean_square(z_hel-zns),...
            root_mean_square(z_helTz-z_hel),  root_mean_square(z_helSz-z_hel),...
            root_mean_square(zns_s2xy-zns), root_mean_square(z_hels2xy-z_hel)];
       
txt1 = '(a) z[$\\omega_+$], RMS = %.2d m';
txt2 = '(b) z[$\\sigma]$ - z[$\\omega_+$], RMS = %6.2f m';
txt3 = '(c) z[$\\omega_s]$ - z[$\\omega_+$], RMS = %6.2f m';
txt4 = '(d) z[$\\omega_{hel}]$ - z[$\\omega_+$], RMS = %6.2f m';
txt5 = '(e) z[$\\omega_{hel\\Theta_z}]$ - z[$\\omega_{hel}$], RMS = %6.2f m';
txt6 = '(f) z[$\\omega_{helS_z}]$ - z[$\\omega_{hel}$], RMS = %6.2f m';
txt7 = '(g) z[$\\omega_{s^2}]$ - z[$\\omega_+$], RMS = %6.2f m';
txt8 = '(h) z[$\\omega_{hel+s^2}]$ - z[$\\omega_{hel}$], RMS = %6.2f m';
title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 6: Streamline
OPTS_AXES = {'Margin', .08, 'Spacing', .08};
figure('Position', [0 0 1800 1000])

x_quiver = repmat(g.XCvec, [1 160]);
y_quiver = repmat(g.YCvec, [360 1]);
u_hel = detail_hel.u;
v_hel = detail_hel.v;
shift = [20 0];

OPTS_FIGS.XSH = 0;
ax1 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, ehel_hel, land_OCCA, OPTS_FIGS);
hold on
quiver(x_quiver, y_quiver, u_hel, v_hel,5, 'k')
hold off
colorbar(ax1)
caxis([-5e-10, 5e-10])
colormap(ax1, bluewhitered), colorbar
txt = '$e^{hel}$ on $\\omega_{hel}$-surface with streamlines';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
ax1.FontSize = 20;

%% Figure 9: The sum of the advection and diffusion
ad_df_ns       = 1e-3*df_ns + ehelns;
ad_df_ns_s     = 1e-3*df_ns_s+ ehelns_s;
ad_df_pt       = 1e-3*df_pt + ehelpt;
ad_df_hel      = 1e-3*df_hel + ehel_hel;
ad_df_helTz    = 1e-3*df_helTz + ehel_helTz;
ad_df_helSz    = 1e-3*df_helSz + ehel_helSz;
ad_df_s2xy     = 1e-3*df_s2xy + ehelns_s2xy;
ad_df_hel_s2xy = 1e-3*df_hel_s2xy + ehel_hel_s2xy;

figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -1e-7*ones(1,8);
OPTS_plot.p_lim = 1e-7*ones(1,8);
OPTS_plot.bwr = ones(1,8);

data = cat(3, ad_df_pt, ad_df_ns, ad_df_ns_s, ad_df_hel,...
                                          ad_df_helTz, ad_df_helSz, ad_df_s2xy, ad_df_hel_s2xy);
data_rms = [root_mean_square(ad_df_pt(:)),    root_mean_square(ad_df_ns(:)),...
            root_mean_square(ad_df_ns_s(:)),  root_mean_square(ad_df_hel(:)),...
            root_mean_square(ad_df_helTz(:)), root_mean_square(ad_df_helSz(:)),...
            root_mean_square(ad_df_s2xy(:)),  root_mean_square(ad_df_hel_s2xy(:))];
quantity_txt = '$e^{hel} + 10^{-3}*D^f$ on ';        
rms_txt = 'RMS  = %.2d m/s';        
txt1 = ['(a) ', quantity_txt,'$\\sigma$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{hel}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{hel\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{helS_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{s^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{hel+s^2}$-surface, ', rms_txt];
title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 10: The scatter of e^hel and s^2 for omega_hel and omega_hel+s^2
figure('Position', [0 0 1200 1000])

OPTS_plot.nbins = 400;

data1 = cat(3, ehelpt, ehelns, ehelns_s, ehel_hel, ehel_helTz, ehel_helSz, ehelns_s2xy, ehel_hel_s2xy);
data1_rms = [ehelpt_rms, ehelns_rms, ehelns_s_rms, ehel_hel_rms, ehel_helTz_rms,...
                                                       ehel_helSz_rms, ehelns_s2xy_rms, ehel_hel_s2xy_rms];

data2 = cat(3, df_pt*1e-3, df_ns*1e-3, df_ns_s*1e-3, df_hel*1e-3,...
                                   df_helTz*1e-3, df_helSz*1e-3, df_s2xy*1e-3, df_hel_s2xy*1e-3);
data2_rms = [s2pt_rms, s2ns_rms, s2ns_s_rms, s2hel_rms, s2helTz_rms,...
                                                       s2helSz_rms, s2ns_s2xy_rms, s2hel_s2xy_rms]; 
txt_x = '$log_{10} (\\ |e^{hel}\\ |)$ [m/s], RMS = %.2d m/s';
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';

txt1 = ['(a) ','$\\sigma$-surface'];
txt2 = ['(b) ','$\\omega_+$-surface'];
txt3 = ['(c) ','$\\omega_s$-surface'];
txt4 = ['(d) ','$\\omega_{hel}$-surface'];
txt5 = ['(e) ','$\\omega_{hel\\Theta_z}$-surface'];
txt6 = ['(f) ','$$\\omega_{helS_z}$-surface'];
txt7 = ['(g) ','$\\omega_{s^2}$-surface'];
txt8 = ['(h) ','$\\omega_{hel+s^2}$-surface'];
title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);

figure_hf = scatter_plotting(data1, data2, data1_rms, data2_rms,txt_x, txt_y, title_text,OPTS_plot, OPTS_FIGS);

