%% salinity contour and velocities on omega+
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])

x_quiver = repmat(g.XCvec, [1 160]);
y_quiver = repmat(g.YCvec, [360 1]);

ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, detailns.u, land_OCCA, OPTS_FIGS);
colorbar(ax3)
% caxis([-1e-7, 1e-7])
hold on
contour(g.XCvec, g.YCvec, sns', 100, 'k','ShowText','on')
quiver(x_quiver, y_quiver, detailns.u_meanT, detailns.v_meanT,5, 'k')
txt = '(a) Salinity Contour and velocity field on $\\omega_+$-surface, mean = %6.2f';
title(sprintf(txt ,nanmean(sns(:))) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);


%% salinity contour and velocities on omega_hel

OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])

ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, s_hel, land_OCCA, OPTS_FIGS);
colorbar(ax3)
% caxis([-1e-7, 1e-7])
hold on
contour(g.XCvec, g.YCvec, s_hel', 100, 'k','ShowText','on')
quiver(x_quiver, y_quiver, detail_hel.u_meanT, detail_hel.v_meanT,5, 'k')
txt = '(a) Salinity Contour and velocity field on $\\omega_{hel}$-surface, mean = %6.2f';
title(sprintf(txt ,nanmean(s_hel(:))) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% salinity contour and velocities on potential density surface

OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])

ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, detailpt.u, land_OCCA, OPTS_FIGS);
colorbar(ax3)
% caxis([-1e-7, 1e-7])
hold on
contour(g.XCvec, g.YCvec, s_sigma', 200, 'k','ShowText','on')
quiver(x_quiver, y_quiver, detailpt.u_meanT, detailpt.v_meanT,5, 'k')
txt = '(a) Salinity Contour and velocity field on $\\omega_{hel}$-surface, mean = %6.2f';
title(sprintf(txt ,nanmean(s_sigma(:))) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% slope error change
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])

ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, 0*s_hel, land_OCCA, OPTS_FIGS);
colormap(ax3, bluewhitered)
quiver(x_quiver, y_quiver, A_X_f(detailns.sx), A_Y_f(detailns.sy),10, 'r')
hold on
quiver(x_quiver, y_quiver, A_X_f(detail_hel.sx), A_Y_f(detail_hel.sy),10, 'k')
title('Slope error change from $\omega_+$ (red) to $\omega_{hel}$ (black)', 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% velocity change
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])

ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, 0*s_hel, land_OCCA, OPTS_FIGS);
colormap(ax3, bluewhitered)
quiver(x_quiver, y_quiver, A_X_f(detailns.u), A_Y_f(detailns.v),5, 'r')
hold on
quiver(x_quiver, y_quiver, A_X_f(detail_hel.u), A_Y_f(detail_hel.v),5, 'k')
title('Velocity change from $\omega_+$ (red) to $\omega_{hel}$ (black)', 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% salinity and temperature contour and gradient on omega+
snsx_grad = D_X_b(sns)./g.DXCvec; % on the U grid
snsy_grad = D_Y_b(sns)./g.DYCsc; % on the V grid

tnsx_grad = D_X_b(tns)./g.DXCvec; % on the U grid
tnsy_grad = D_Y_b(tns)./g.DYCsc; % on the V grid

OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])

ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, 0*s_hel, land_OCCA, OPTS_FIGS);
colormap(ax3, bluewhitered)
hold on
contour(g.XCvec, g.YCvec, sns', 100, 'k')
contour(g.XCvec, g.YCvec, tns', 100, 'r')
quiver(x_quiver', y_quiver', A_X_f(snsx_grad)', A_Y_f(snsy_grad)',10, 'k')
quiver(x_quiver', y_quiver', A_X_f(tnsx_grad)', A_Y_f(tnsy_grad)',10, 'r')
title('Salinity contour (with $\nabla_a S$) (red) and potential temperature contour (with $\nabla_a \Theta$) (black) on $\omega_+$', 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% salinity and temperature contour and gradient on omega_hel
thelx_grad = D_X_b(t_hel)./g.DXCvec; % on the U grid
thely_grad = D_Y_b(t_hel)./g.DYCsc; % on the V grid

shelx_grad = D_X_b(s_hel)./g.DXCvec; % on the U grid
shely_grad = D_Y_b(s_hel)./g.DYCsc; % on the V grid


OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])
ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, 0*s_hel, land_OCCA, OPTS_FIGS);
colormap(ax3, bluewhitered)
hold on
contour(g.XCvec, g.YCvec, s_hel', 100, 'k')
contour(g.XCvec, g.YCvec, t_hel', 100, 'r')
quiver(x_quiver', y_quiver', A_X_f(shelx_grad)', A_Y_f(shely_grad)',10, 'k')
quiver(x_quiver', y_quiver', A_X_f(thelx_grad)', A_Y_f(thely_grad)',10, 'r')
title('Salinity contour (with $\nabla_a S$) (red) and potential temperature contour (with $\nabla_a \Theta$) (black) on $\omega_{hel}$', 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% salinity gradient change from omega+ to omega_hel
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])
ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, 0*s_hel, land_OCCA, OPTS_FIGS);
colormap(ax3, bluewhitered)
hold on
quiver(x_quiver', y_quiver', A_X_f(snsx_grad)', A_Y_f(snsy_grad)',10, 'r')
quiver(x_quiver', y_quiver', A_X_f(shelx_grad)', A_Y_f(shely_grad)',10, 'k')
title('$\nabla_a S$ from $\omega_+$ (red) to $\omega_{hel}$ (black)' , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% salinity gradient change from omega+ to omega_hel (merge)
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])
ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, 0*s_hel, land_OCCA, OPTS_FIGS);
colormap(ax3, bluewhitered)
hold on
quiver(x_quiver', y_quiver', A_X_f(snsx_grad)' - A_X_f(shelx_grad)', A_Y_f(snsy_grad)' - A_Y_f(shely_grad)',10, 'r')
title('$\nabla_a S$ from $\omega_+$ (red) to $\omega_{hel}$ (black)' , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% T gradient change from omega+ to omega_hel
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])
ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, 0*s_hel, land_OCCA, OPTS_FIGS);
colormap(ax3, bluewhitered)
hold on
quiver(x_quiver', y_quiver', A_X_f(snsx_grad)', A_Y_f(snsy_grad)',10, 'r')
quiver(x_quiver', y_quiver', A_X_f(thelx_grad)', A_Y_f(thely_grad)',10, 'k')
title('$\nabla_a \Theta$ from $\omega_+$ (red) to $\omega_{hel}$ (black)' , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% T gradient change from omega+ to omega_hel (merge)
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])
ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, 0*s_hel, land_OCCA, OPTS_FIGS);
colormap(ax3, bluewhitered)
hold on
quiver(x_quiver', y_quiver', A_X_f(snsx_grad)' - A_X_f(thelx_grad)', A_Y_f(snsy_grad)' - A_Y_f(thely_grad)',10, 'r')
title('$\nabla_a \Theta$ from $\omega_+$ (red) to $\omega_{hel}$ (black)' , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%% rho_s and rho_\Theta on the omega+
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])

[rs_ns, rt_ns] = densjmd95_bsq_second_derivs(sns, tns, zns);
[rs_hel, rt_hel] = densjmd95_bsq_second_derivs(s_hel, t_hel, z_hel);
rs_sx_ns = rs_ns.* A_X_f(snsx_grad);

ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, rs_sx_ns, land_OCCA, OPTS_FIGS);
colorbar
txt = '$\\rho_S S_x$ on $\\omega_+$-surface, mean = %.2d';
title(sprintf(txt ,root_mean_square(rs_sx_ns(:))) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%%
[ex_ns,ey_ns] = ntp_errors(sns, tns, zns, g.DXCvec, g.DYCsc, 1, 0, [1;1]);

OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])
ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, ex_ns, land_OCCA, OPTS_FIGS);
colorbar
txt = '$\\epsilon_x$ on $\\omega_+$-surface, mean = %.2d';
title(sprintf(txt ,root_mean_square(ex_ns(:))) , 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

%%
OPTS_AXES = {'Margin', .05, 'Spacing', .05};
figure('Position', [0 0 1800 1000])

ax3 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, rt_ns, land_OCCA, OPTS_FIGS);
colorbar
title('$\rho_{\Theta}$ on the $\omega_+$-surface', 'fontsize',10,'Interpreter','latex');
set(gca,'fontsize', 18);

