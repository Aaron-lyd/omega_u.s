%% load data
ni = 2^6+1; nj = ni; nk = 20; wall = 1; helicity = 5; U_M = 0;
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
U_smooth = smooth2a(U_rand, 20,10,[],[1,1]);

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
OPTS_DC.FIGS_SHOW = 0;

% velocity
v_mode =2;
if v_mode ==1 % purely rotational random
    u_field = squeeze(U_smooth(1,:,:))';
    U_f = U_curl;
    V_f = V_curl;
    OPTS_DC.WRAP = [1;1];
elseif v_mode ==2 % purely irrotational random
    u_field = squeeze(U_smooth(1,:,:))';
    U_f = U_grad;
    V_f = V_grad;
    OPTS_DC.WRAP = [1;1];
elseif v_mode ==3 % irrotational Gaussian (on checkerboard)
    u_field = squeeze(phii(1,:,:))';
    U_f = 10*U_M;
    V_f = 10*V_M;
    OPTS_DC.WRAP = [1;1];
end

OPTS_DC.PCG = 0;
OPTS_DC.ITER = 50;
OPTS_DC.SHEAR = 1;
OPTS_DC.STRAT = 1;
OPTS_DC.CHOLESKY = 0;
OPTS_DC.H_SIM = 0;

OPTS_DC.LM = 0;
% OPTS_DC.TK = 1e-5;
OPTS_DC.TK = 1e-8;

OPTS_DC.MODE = 1; %omega_u.s surface

[z_hel_dc,s_hel_dc,t_hel_dc, ~,d_hel_dc] = omega_hel_surface(S, T, P, U_f, V_f, zns_dc, OPTS_DC);

%% omega_s
OPTS_DC.MODE = 2;
OPTS_DC.STRAT = 1;
OPTS.SHEAR = 1;
% OPTS.LM = 1e-3;
OPTS_DC.TK = 0;
OPTS_DC.H_SIM = 0;
% OPTS.LM = 0;
% OPTS.TK = 1e-8;
OPTS_DC.TOL_LSQR_REL = 1e-6;
OPTS_DC.CHOLESKY = 0;

OPTS_DC.ITER = 100;

[zns_s_dc,sns_s_dc,tns_s_dc, ~, d_s_dc] = omega_hel_surface(S, T, P, U_f, V_f, zns_dc, OPTS_DC);

%% omega_s2xy
OPTS.MODE = 6;
OPTS.STRAT = 1;
OPTS.SHEAR = 1;
OPTS.ITER = 20;
OPTS.H_SIM = 0;
OPTS.CHOLESKY = 0;


OPTS.LM = 0;
% OPTS.TK = 1e-9;
OPTS.TK = 1e-9;

% OPTS.TK = 0;
% OPTS.LM = 1;

[zns_s2xy_dc,sns_s2xy_dc,tns_s2xy_dc, ~, d_s2xy_dc] = omega_hel_surface(S, T, P, U_f, V_f, zns_dc, OPTS_DC);

%% omega_hel_s2xy 
OPTS.ITER = 100;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;
OPTS.CHOLESKY = 0;

OPTS.LM = 0;
OPTS.TK = 1e-8;
OPTS.H_SIM = 0;

% OPTS.TK = 0;
% OPTS.LM = 1;

% weight of ehelu and ehelv
OPTS.s2xyW = 1;
OPTS.FIGS_SHOW = 0;

OPTS.MODE = 8;

[z_hels2xy_dc,s_hels2xy_dc,t_hels2xy_dc, ~,d_hels2xy_dc] = omega_hel_surface(S, T, P, U_f, V_f, zns_dc, OPTS_DC);

%% omega_hel_Tz
OPTS.data_cube=0;
OPTS.ITER = 100;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;
OPTS.CHOLESKY = 0;

OPTS.LM = 0;
OPTS.TK = 1e-9;

OPTS.MODE = 9;

[z_helTz_dc,s_helTz_dc,t_helTz_dc, ~,d_helTz_dc] = omega_hel_surface(S, T, P, U_f, V_f, zns_dc, OPTS_DC);

%% omega_hel_Sz
OPTS.data_cube=0;
OPTS.ITER = 100;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;
OPTS.CHOLESKY = 0;

OPTS.LM = 0;
OPTS.TK = 1e-10;

OPTS.MODE = 10;

[z_helSz_dc,s_helSz_dc,t_helSz_dc, ~,d_helSz_dc] = omega_hel_surface(S, T, P, U_f, V_f, zns_dc, OPTS_DC);

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
[s_ns_dc, ~, s_ns_dc_rms, ~, ehelns_dc, ehelns_dc_rms, ~, detailns_dc] =...
    s_ehel_flux(sns_dc,tns_dc, zns_dc, S, T, P, U_f, V_f, [1;1], OPT_DC);

% s, ehel and flux for omega_hel
OPT_DC.MODES = 0; % using NTP to calculate slope errors
OPT_DC.MODEHEL = 1; % using divergence method to calculate e^hel
[shel_dc, ~, s_hel_dc_rms, ~, ehel_hel_dc, ehel_hel_dc_rms, ~, detail_hel_dc] =...
    s_ehel_flux(s_hel_dc,t_hel_dc, z_hel_dc, S, T, P, U_f, V_f, [1;1], OPT_DC);

% cos theta
costheta_ns = ehelns_dc./(abs(s_ns_dc).*abs(detailns_dc.u_t));
costheta_hel = ehel_hel_dc./(abs(shel_dc).*abs(detail_hel_dc.u_t));

% Helicity
A_vor = g.RAZvec(80); % the area of the vorticity cell
% A_vor_3D = permute(repmat(A_vor, [360 1 50]), [3 1 2]);

H_3D = hel(S, T, P, A_vor);

H_hel = interpfn(P, H_3D, z_hel_dc);
H_ns = interpfn(P, H_3D, zns_dc);

rms(H_hel(:))
rms(H_ns(:))


%% Figure 1:e^hel maps
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1400 1000])

ax1 = subaxis(3,2,1, OPTS_AXES{:});
imagesc(u_field);
hold on
contour(u_field, 40, 'k')
quiver(A_X(squeeze(U_f(1,:,:)))', A_Y(squeeze(V_f(1,:,:)))', 1, 'k')
title('Velocity field', 'fontsize',16,'Interpreter','latex');

ax2 = subaxis(3,2,2, OPTS_AXES{:});
pcolor(ehel_hel_dc)
shading flat
%caxis([-5e-13, 5e-13])
colormap(ax2, bluewhitered), colorbar
X = (1:ni)';
Y = (1:nj);
pcol_logsigned(ax2, X, Y, ehel_hel_dc, -19, -14);
txt = '$e^{u}$ on $\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface, RMS = %.2d m/s';
title(sprintf(txt ,ehel_hel_dc_rms) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

ax3 = subaxis(3,2,3, OPTS_AXES{:});
pcolor(log10(shel_dc))
shading flat
%caxis([-5e-13, 5e-13])
% colormap(ax3, bluewhitered), 
colorbar
txt = '$log_{10} (|\\bf{s}|)$ on $\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface, RMS = %.2d m/s';
title(sprintf(txt ,s_hel_dc_rms) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

ax4 = subaxis(3,2,4, OPTS_AXES{:});
pcolor(costheta_hel)
shading flat
caxis([-1, 1])
cm = parula(64);
cm = [cm(1:end-1,:); flipud(cm)];
colormap(ax4,cm);
colorbar(ax4)
txt = '$\\mathbf{u} \\cdot \\mathbf{s}/|\\mathbf{u}||\\mathbf{s}|$ on $\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

ax5 = subaxis(3,2,5, OPTS_AXES{:});
pcolor(H_hel)
shading flat
%caxis([-5e-13, 5e-13])
colormap(ax5, bluewhitered), colorbar
pcol_logsigned(ax5, X, Y, H_hel, -25, -19);
txt = 'Helicity on $\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% Figure 2: maps for and omega+
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1600 1000])

ax1 = subaxis(2,2,1, OPTS_AXES{:});
pcolor(ehelns_dc)
shading flat
% caxis([-5e-14, 5e-14]) 
colormap(ax1, bluewhitered), colorbar
pcol_logsigned(ax1, X, Y, ehelns_dc, -15, -8);

txt = '(a) $e^{u}$ on $\\omega_+$-surface, RMS = %.2d m/s';
title(sprintf(txt ,ehelns_dc_rms) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

ax2 = subaxis(2,2,2, OPTS_AXES{:});
pcolor(detail_hel_dc.sx)
shading flat
caxis([-5e-8, 5e-8])
% colormap(ax3, bluewhitered), 
colorbar
txt = '$log_{10} (|\\bf{s}|)$ on $\\omega_+$-surface, RMS = %.2d m/s';
title(sprintf(txt ,s_ns_dc_rms) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

ax3 = subaxis(2,2,3, OPTS_AXES{:});
pcolor(costheta_ns)
shading flat
caxis([-1, 1])
cm = parula(64);
cm = [cm(1:end-1,:); flipud(cm)];
colormap(ax3,cm);
colorbar(ax3)
txt = '$\\mathbf{u} \\cdot \\mathbf{s}/|\\mathbf{u}||\\mathbf{s}|$ on $\\omega_+$-surface';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

ax4 = subaxis(2,2,4, OPTS_AXES{:});
pcolor(H_ns)
shading flat
%caxis([-5e-13, 5e-13])
colormap(ax4, bluewhitered), colorbar
pcol_logsigned(ax4, X, Y, H_ns, -25, -19);
txt = 'Helicity on $\\omega_+$-surface';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);


%%
figure
pcolor(detail_hel_dc.sx)
shading flat
caxis([-5e-8, 5e-8])
colorbar
txt = '$s_x$ on $\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

