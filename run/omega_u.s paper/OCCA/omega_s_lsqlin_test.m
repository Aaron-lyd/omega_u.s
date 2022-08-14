%% omega_s
OPTS.MODE = 2;
OPTS.STRAT = 1;
OPTS.CHOLESKY = 0;
OPTS.TK = 0;
OPTS.LM = 0;
OPTS.H_SIM = 0;
OPTS.KMJ = 0;
OPTS.FIGS_SHOW = 1;
OPTS.ITER = 50;

[zns_s,sns_s,tns_s, ~, d_s] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%
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

% s, ehel and flux for omega_s
OPT.MODES = 0; % using NTP to calculate slope errors
OPT.MODEHEL = 1; % using divergence method to calculate e^hel
[s_ns_s, s2ns_s, sns_s_rms, s2ns_s_rms, ehelns_s, ehelns_s_rms, fluxns_s, detailns_s] =...
                             s_ehel_flux(sns_s,tns_s, zns_s, SB, TB, ZB, U, V, [1;1], OPT);

%% slope error
OPTS_AXES = {'Margin', .04, 'Spacing', .03};
figure('Position', [0 0 1800 1000])

ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec,  (log10(s_ns_s)), land_OCCA, OPTS_FIGS);
colorbar(ax3)
txt = '$log_{10} \\bf(s)$  on $\\omega_s$-surface, RMS = %.2d';
title(sprintf(txt ,sns_s_rms) , 'fontsize',10,'Interpreter','latex');
caxis([-8, -3])
set(gca,'fontsize', 15);
ax3.XTickLabel = [];

%% depth difference
OPTS_AXES = {'Margin', .04, 'Spacing', .03};
figure('Position', [0 0 1800 1000])
ax2 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax2, g.XCvec, g.YCvec,  (zns_s-zns), land_OCCA, OPTS_FIGS);
caxis([-10,10])
colormap(ax2,bluewhitered), colorbar
txt = 'z[$\\omega_s]$ - z[$\\omega_+$], RMS = %6.2f m';
title(sprintf(txt ,root_mean_square(zns_s-zns)) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 15);

[~,i0] = min(abs(g.XCvec - OPTS.x0));
[~,j0] = min(abs(g.YCvec - OPTS.y0));

pinning_z_diff = zns_s(i0,j0) -  zns(i0, j0);
