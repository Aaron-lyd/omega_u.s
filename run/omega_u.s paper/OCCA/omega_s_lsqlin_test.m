%% omega_s
OPTS.MODE = 11;
OPTS.STRAT = 1;
OPTS.SHEAR = 1;
% OPTS.LM = 1e-3;
OPTS.TK = 0;
OPTS.H_SIM = 0;
% OPTS.LM = 0;
% OPTS.TK = 1e-8;
OPTS.TOL_LSQR_REL = 1e-6;

OPTS.ITER = 100;

[zns_s,sns_s,tns_s, ~, d_s] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%%
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

%%
OPTS_AXES = {'Margin', .04, 'Spacing', .03};
figure('Position', [0 0 1200 1000])

ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec,  (log10(s_ns_s)), land_OCCA, OPTS_FIGS);
colorbar(ax3)
txt = '$log_{10} \\bf(s)$  on $\\omega_s$-surface, RMS = %.2d';
title(sprintf(txt ,sns_s_rms) , 'fontsize',10,'Interpreter','latex');
caxis([-8, -3])
set(gca,'fontsize', 15);
ax3.XTickLabel = [];