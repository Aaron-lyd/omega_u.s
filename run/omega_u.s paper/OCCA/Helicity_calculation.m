%%
A_vor = g.RAZvec; % the area of the vorticity cell
A_vor_3D = permute(repmat(A_vor, [360 1 50]), [3 1 2]);

H_3D = hel(SB, TB, ZB, A_vor_3D);

H_hel = interpfn(ZB, H_3D, z_hel);
H_ns = interpfn(ZB, H_3D, zns);

%%

OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1200 1000])


ax2 = subaxis(2,1,1, OPTS_AXES{:});
hf = fig_map(ax2, g.XCvec, g.YCvec, H_ns, land_OCCA, OPTS_FIGS);
colorbar(ax2)
caxis([-1e-16, 1e-16])
colormap(ax2, bluewhitered), colorbar
txt = '(a) Helicity on $\\omega_{+}$-surface, RMS = %.4d';
title(sprintf(txt ,root_mean_square(H_ns)) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);
ax2.XTickLabel = [];

ax2 = subaxis(2,1,2, OPTS_AXES{:});
hf = fig_map(ax2, g.XCvec, g.YCvec, H_hel, land_OCCA, OPTS_FIGS);
colorbar(ax2)
caxis([-1e-16, 1e-16])
colormap(ax2, bluewhitered), colorbar
txt = '(b) Helicity on $\\omega_{hel}$-surface, RMS = %.4d';
title(sprintf(txt ,root_mean_square(H_hel)) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);
ax2.XTickLabel = [];

%% Hel_diff
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])


ax2 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax2, g.XCvec, g.YCvec, H_ns - H_hel, land_OCCA, OPTS_FIGS);
colorbar(ax2)
caxis([-1e-18, 1e-18])
colormap(ax2, bluewhitered), colorbar
txt = 'H[$\\omega_{+}$-surface] - H[$\\omega_{hel}$-surface], RMS = %.4d';
title(sprintf(txt ,root_mean_square(H_ns - H_hel)) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);
ax2.XTickLabel = [];
