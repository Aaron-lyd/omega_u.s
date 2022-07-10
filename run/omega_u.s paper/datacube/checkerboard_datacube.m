%% load data
ni = 2^6+1; nj = ni; nk = 20; wall = 1; helicity = 0.2; U_M = 0;
[S, T, P,U_M,V_M, ~, phii] = datacube(ni,nj,nk,helicity,wall,U_M,1,0);

% interpolate salinity and temperature on the neutral density surface
% interpfn = @ppc_linterp;
interpfn = @ppc_pchip;
SppZ = interpfn(P, S);
TppZ = interpfn(P, T);

lead1 = @(x) reshape(x, [1 size(x)]);

% reference cast
i0 = floor(ni/2)+1;
j0 = floor(ni/2)+1;
I0 = sub2ind(size(squeeze(S(1,:,:))), i0, j0);

% Start with a potential density surface, crudely calculated
z_ref = 500;
z0 = 1000;

SIGMA = eos(S, T, z_ref); % the potential density of the world 3D data reference to z_ref
sigma0 = interpfn(P, SIGMA(:,i0,j0), z0); % The isovalue of the potential density we want
z_pt = interpfn(SIGMA, P, sigma0); % the potential density surface's depth of the potential density surface reference to z_ref with iso value sigma0
s_pt = interpfn(P, S, z_pt);
t_pt = interpfn(P, T, z_pt);

%% Codegen for ppc_val* and ntp_slope
[nk, ni, nj] = size(S);
t_Z = coder.typeof(0, [nk, 1], [0, 0]);
t_Sppc = coder.typeof(0, [4, nk-1, ni, nj], [1, 0, 0, 0]);
t_s = coder.typeof(0, [ni, nj], [1, 1]);
codegen('ppc_val', '-report', '-args', {t_Z, t_Sppc, t_s, 0})
codegen('ppc_val2', '-report', '-args', {t_Z, t_Sppc, t_Sppc, t_s, 0})
ni_ = max(4096, ni);
nj_ = max(4096, nj);
nk_ = max(128, nk);
ntp_slope_codegen(nk_, ni_, nj_, isvector(P));
ntp_slope_codegen(nk, ni, nj, isvector(P));

%% define some functions
im1 = @(F) circshift(F, [0 +1 0]);
jm1 = @(F) circshift(F, [0 0 +1]);
ip1 = @(F) circshift(F, [0 -1 0]);
jp1 = @(F) circshift(F, [0 0 -1]);

% avergece and backward difference
A_X = @(F) (F + im1(F)) / 2;
D_X_f = @(F)  ip1(F) - F;
A_Y = @(F) (F + jm1(F)) / 2;
D_Y_f = @(F)  jp1(F) - F;

A_X_f = @(F) (F + ip1(F)) / 2;
A_Y_f = @(F) (F + jp1(F)) / 2;

D_X_b = @(F)  F - im1(F);
D_Y_b = @(F)  F - jm1(F);

%% fake velocity
% build a 2D smooth velocity field
U_rand = randn(ni, nj);
U_smooth = smooth2a(U_rand, 10,3,[],[1,1]);

% make it to be 3D
U_smooth = permute(repmat(U_smooth, [1 1 nk]), [3 1 2]); % on vorticity/tracer grid

% dx and dy
dx_3D = repmat(g.DXCvec(80), [nk ni nj]);
dy_3D = repmat(g.DYCsc, [nk ni nj]);

% rotational, divergence free velocity field, checkerboard everywhere
U_curl = -D_Y_f(U_smooth); % on the U grid
V_curl = D_X_f(U_smooth);  % on the V grid

% purely divergence vector velocity field, no checkerboard
U_grad = D_X_b(U_smooth); % on the U grid
V_grad = D_Y_b(U_smooth); % on the V grid

%% velocity quiver for purely rotational velocotiy field
OPTS_AXES = {'Margin', .05, 'Spacing', .05};

figure('Position', [0 0 1800 1000])
imagesc(squeeze(U_smooth(1,:,:))');
hold on
contour(squeeze(U_smooth(1,:,:))', 40, 'k')
quiver(A_X_f(squeeze(U_curl(1,:,:)))', A_Y_f(squeeze(V_curl(1,:,:)))', 1, 'k')
title('Rotational velocity field', 'fontsize',16,'Interpreter','latex');

%% velocity quiver for purely irrotational velocotiy field
OPTS_AXES = {'Margin', .05, 'Spacing', .05};

figure('Position', [0 0 1800 1000])
imagesc(squeeze(U_smooth(1,:,:))');
hold on
contour(squeeze(U_smooth(1,:,:))', 40, 'k')
quiver(A_X_f(squeeze(U_grad(1,:,:)))', A_Y_f(squeeze(V_grad(1,:,:)))', 1, 'k')
title('Irrotational velocity field', 'fontsize',16,'Interpreter','latex');

%% velocity quiver for Gaussian velocotiy field
OPTS_AXES = {'Margin', .05, 'Spacing', .05};

figure('Position', [0 0 1800 1000])
imagesc(squeeze(phii(1,:,:))');
colorbar
hold on
contour(squeeze(phii(1,:,:))', 40, 'k')
quiver(A_X_f(squeeze(U_M(1,:,:)))', A_Y_f(squeeze(V_M(1,:,:)))', 1, 'k')
title('Gaussian velocity field', 'fontsize',16,'Interpreter','latex');

%% omega_1
OPTS_DC = [];
OPTS_DC.INTERPFN = interpfn;

[zns_dc, sns_dc, tns_dc, dns_dc] = omega_surface(S, T, P, z_pt, I0,[1;1], OPTS_DC); %in x-y direction

%% Grid information
OPTS_DC.data_cube = 1;
OPTS_DC.nx = ni; 
OPTS_DC.ny = nj;
OPTS_DC.i0 = i0;
OPTS_DC.j0 = j0;

OPTS_DC.DXCvec = g.DXCvec(80); %dx
OPTS_DC.DXGvec = g.DXGvec(80);
OPTS_DC.DYCsc = g.DYCsc; %dy
OPTS_DC.DYGsc = g.DYGsc;
OPTS_DC.RACvec = g.RACvec(80);% Vertical area of the tracer cells [m^2]

% OPTS set up
OPTS_DC.DAMP = 1;
OPTS_DC.INTERPFN = @ppc_pchip;
% OPTS.WRAP = [1;1];
OPTS_DC.FIGS_SHOW = 1;

% omega_hel
v_mode =2;
if v_mode ==1 % checkerboard
    U_f = U_curl;
    V_f = V_curl;
    OPTS_DC.WRAP = [1;1];
elseif v_mode ==2 % no checkerboard
    U_f = U_grad;
    V_f = V_grad;
    OPTS_DC.WRAP = [1;1];
elseif v_mode ==3
    U_f = U_M;
    V_f = V_M;
    OPTS_DC.WRAP = [1;1];
end

OPTS_DC.ITER = 1;
OPTS_DC.SHEAR = 0;
OPTS_DC.STRAT = 0;

OPTS_DC.LM = 0;
%OPTS_DC.TK = 1e-5;
OPTS_DC.TK = 0;

OPTS_DC.MODE = 1;

[z_hel_dc,s_hel_dc,t_hel_dc, ~,d_hel_dc] = omega_hel_surface(S, T, P, U_f, V_f, zns_dc, OPTS_DC);

%% s, e^hel and flux
% prepare for the grid information
OPT_DC.DXCvec = g.DXCvec(80); % dx
OPT_DC.DYCsc = g.DYCsc; % dy
OPT_DC.DYGsc = g.DYGsc;
OPT_DC.DXGvec = g.DXGvec(80);
OPT_DC.RAWvec = g.RAWvec(80); % area of the grid in u grid
OPT_DC.RASvec = g.RASvec(80); % area of the grid in v grid
OPT_DC.RACvec = g.RACvec(80); % area of the grid in tracer grid
OPT_DC.nx = ni;
OPT_DC.ny = nj;
OPT_DC.tolz = 1e-8; % tolerance of slope error

% s, ehel and flux for omega+
OPT_DC.MODES = 0; % using epsilon to calculate slope errors
OPT_DC.MODEHEL = 1; % Using divergence method to calculate e^hel
[~, ~, ~, ~, ehelns_dc, ehelns_dc_rms, ~, detailns] =...
    s_ehel_flux(sns_dc,tns_dc, zns_dc, S, T, P, U_f, V_f, [1;1], OPT_DC);

% s, ehel and flux for omega_hel
OPT_DC.MODES = 0; % using NTP to calculate slope errors
OPT_DC.MODEHEL = 1; % using divergence method to calculate e^hel
[~, ~, ~, ~, ehel_hel_dc, ehel_hel_dc_rms, ~, ~] =...
    s_ehel_flux(s_hel_dc,t_hel_dc, z_hel_dc, S, T, P, U_f, V_f, [1;1], OPT_DC);

%% Figure 1:e^hel maps
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1000 1000])

ax1 = subaxis(2,1,1, OPTS_AXES{:});
imagesc(squeeze(U_smooth(1,:,:))');
hold on
contour(squeeze(U_smooth(1,:,:))', 40, 'k')
quiver(A_X(squeeze(U_f(1,:,:)))', A_Y(squeeze(V_f(1,:,:)))', 1, 'k')
title('Velocity field', 'fontsize',16,'Interpreter','latex');

ax2 = subaxis(2,1,2, OPTS_AXES{:});
pcolor(ehel_hel_dc)
shading flat
% caxis([-5e-13, 5e-13])
colormap(ax2, bluewhitered), colorbar
txt = '$e^{hel}$ on $\\omega_{hel}$-surface, RMS = %.2d m/s';
title(sprintf(txt ,ehel_hel_dc_rms) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% Figure 1:e^hel maps
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1000 1000])

ax1 = subaxis(2,1,1, OPTS_AXES{:});
pcolor(ehelns_dc)
shading flat
% caxis([-5e-14, 5e-14]) 
colormap(ax1, bluewhitered), colorbar
txt = '(a) $e^{hel}$ on $\\omega_+$-surface, RMS = %.2d m/s';
title(sprintf(txt ,ehelns_dc_rms) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

ax2 = subaxis(2,1,2, OPTS_AXES{:});
pcolor(ehel_hel_dc)
shading flat
% caxis([-5e-14, 5e-14])
colormap(ax2, bluewhitered), colorbar
txt = '(b) $e^{hel}$ on $\\omega_{hel}$-surface, RMS = %.2d m/s';
title(sprintf(txt ,ehel_hel_dc_rms) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);




