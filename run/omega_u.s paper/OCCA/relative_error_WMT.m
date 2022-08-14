%%
OPTS_FIGS.LANDCOL = [1 1 1]*0; % Black
OPTS_FIGS.NANCOL = [1 1 1]* .75; % Grey
OPTS_FIGS.LATLIM = [-80, 72]; % The z(180,0) = 2000 surfaces exclude the Arctic.
OPTS_AXES = {'Margin', .08, 'Spacing', .08};
figure('Position', [0 0 1800 1000])
land_OCCA = squeeze(isnan(SB(1,:,:)));

ax1 = subaxis(2,2,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec,  Tz_ns, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = '$log_{10} (s)$  on $\\omega$-surface (OCCA), RMS = %.2d';
title(sprintf(txt ,root_mean_square(Tz_ns)) , 'fontsize',10,'Interpreter','latex');
%caxis([-8, -3])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%%

dzx_us = (z_hel - im1(z_hel))./g.DXCvec;
dzy_us = (z_hel - jm1(z_hel))./g.DYCsc;

dzx_us2 = dzx_us .* dzx_us .* (g.DXCvec .* g.DYGsc);
dzy_us2 = dzy_us .* dzy_us .* (g.DXGvec .* g.DYCsc);
dz_us2 = (dzx_us2 + ip1(dzx_us2) + dzy_us2 + jp1(dzy_us2)) ./ (2 * g.RACvec);
dz_us = sqrt(dz_us2);

figure('Position', [0 0 1800 1000])
ax1 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec,  log10(dz_us), land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = '$log_{10}(|\\nabla_a z \\cdot \\nabla_a z|)$  on $\\omega_{u.s}$-surface (OCCA), RMS = %.2d';
title(sprintf(txt ,root_mean_square(dz_us)) , 'fontsize',10,'Interpreter','latex');
%caxis([-8, -3])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%%

sTz = Tz_hel .* ss_hel;
grad_hel_T = dotproduct_div(t_hel, g.DXCvec, g.DYCsc, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

r_hel = sTz./ grad_hel_T;

figure('Position', [0 0 1800 1000])
ax1 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, r_hel, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = '$s\\Theta_z/\\nabla_a \\Theta$  on $\\omega_{u.s}$-surface (OCCA), RMS = %.2d';
title(sprintf(txt ,root_mean_square(r_hel)) , 'fontsize',10,'Interpreter','latex');
%caxis([-8, -3])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%%

Pppc = ppc_pchip(ZB, pB);

p_hel = ppc_pchip(ZB, pB, z_hel);

Pz_hel = ppc_val_mex(ZB, Pppc, p_hel, 1);

sPz = Pz_hel .* ss_hel;
grad_hel_P = dotproduct_div(p_hel, g.DXCvec, g.DYCsc, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

r_hel = sPz./ grad_hel_P;

figure('Position', [0 0 1800 1000])
ax1 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, r_hel, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = '$s p_z/\\nabla_a p$  on $\\omega_{u.s}$-surface (OCCA), RMS = %.2d';
title(sprintf(txt ,root_mean_square(r_hel)) , 'fontsize',10,'Interpreter','latex');
%caxis([-8, -3])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%% TB relative error on omega_u.s

grad_p_hel_x = (p_hel - im1(p_hel))./g.DXCvec;
grad_p_hel_y = (p_hel - jm1(p_hel))./g.DYCsc;

grad_T_hel_x = (t_hel - im1(t_hel))./g.DXCvec;
grad_T_hel_y = (t_hel - jm1(t_hel))./g.DYCsc;

s2pzTz = Pz_hel.* Tz_hel .* s2hel;
pzgradTs_x = Pz_hel.*grad_T_hel_x + Tz_hel .* grad_p_hel_x;
pzgradTs_y = Pz_hel.*grad_T_hel_y + Tz_hel .* grad_p_hel_y;

pzgradTs = div_dot(detail_hel.sx, detail_hel.sy,pzgradTs_x ,pzgradTs_y, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

grad_p_grad_T = div_dot(grad_p_hel_x, grad_p_hel_y, grad_T_hel_x, grad_T_hel_y, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

rel_hel = (s2pzTz + pzgradTs)./grad_p_grad_T;

figure('Position', [0 0 1800 1000])
ax1 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, rel_hel, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = 'relative error  on $\\omega_{u.s}$-surface (OCCA), RMS = %.2f';
title(sprintf(txt ,nanrms(rel_hel(:))) , 'fontsize',10,'Interpreter','latex');
caxis([0, 0.1])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%% TB relative error on omega_s

p_ns_s = ppc_pchip(ZB, pB, zns_s);

Pz_ns_s = ppc_val_mex(ZB, Pppc, zns_s, 1);

grad_p_s_x = (p_ns_s - im1(p_ns_s))./g.DXCvec;
grad_p_s_y = (p_ns_s - jm1(p_ns_s))./g.DYCsc;

grad_T_s_x = (tns_s - im1(tns_s))./g.DXCvec;
grad_T_s_y = (tns_s - jm1(tns_s))./g.DYCsc;

s2pzTz = Tz_ns_s.* Pz_ns_s .* s2ns_s;
pzgradTs_x = Pz_ns_s.*grad_T_s_x + Tz_ns_s .* grad_p_s_x;
pzgradTs_y = Pz_ns_s.*grad_T_s_y + Tz_ns_s .* grad_p_s_y;

pzgradTs = div_dot(detailns_s.sx, detailns_s.sy,pzgradTs_x ,pzgradTs_y, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

grad_p_grad_T = div_dot(grad_p_s_x, grad_p_s_y, grad_T_s_x, grad_T_s_y, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

rel_hel = (s2pzTz + pzgradTs)./grad_p_grad_T;

figure('Position', [0 0 1800 1000])
ax1 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, rel_hel, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = 'Thermobaric relative error on $\\omega_{s}$-surface (OCCA), RMS = %.2f';
title(sprintf(txt ,nanrms(rel_hel(:))) , 'fontsize',10,'Interpreter','latex');
caxis([0, 0.1])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%% CB relative error on omega_s

p_ns_s = ppc_pchip(ZB, pB, zns_s);

Pz_ns_s = ppc_val_mex(ZB, Pppc, p_ns_s, 1);

grad_p_s_x = (p_ns_s - im1(p_ns_s))./g.DXCvec;
grad_p_s_y = (p_ns_s - jm1(p_ns_s))./g.DYCsc;

grad_T_s_x = (tns_s - im1(tns_s))./g.DXCvec;
grad_T_s_y = (tns_s - jm1(tns_s))./g.DYCsc;

s2TzTz = Tz_ns_s.* Tz_ns_s .* s2ns_s;
TzgradTs_x = Tz_ns_s .* grad_T_s_x;
TzgradTs_y = Tz_ns_s .* grad_T_s_y;

TzgradTs = div_dot(detailns_s.sx, detailns_s.sy,TzgradTs_x ,TzgradTs_y, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

grad_T_grad_T = div_dot(grad_T_s_x, grad_T_s_y, grad_T_s_x, grad_T_s_y, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

rel_cb_s = (s2TzTz + 2*TzgradTs)./grad_T_grad_T;

figure('Position', [0 0 1800 1000])
ax1 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, rel_cb_s, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = 'Cabeling relative error on $\\omega_{s}$-surface (OCCA), RMS = %.2f';
title(sprintf(txt ,nanrms(rel_cb_s(:))) , 'fontsize',10,'Interpreter','latex');
caxis([0, 0.1])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%%
[TB_rel_s, CB_rel_s] = rel_TB_CB(tns_s, zns_s,detailns_s.sx, detailns_s.sy, s2ns_s, ZB, pB, TB, g.DXCvec, g.DYCsc, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

figure('Position', [0 0 1200 1000])
ax1 = subaxis(2,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, CB_rel_s, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = 'Cabeling relative error on $\\omega_{s}$-surface (OCCA), RMS = %.2f';
title(sprintf(txt ,nanrms(CB_rel_s(:))) , 'fontsize',10,'Interpreter','latex');
caxis([0, 0.1])
ax1.FontSize = 15;
ax1.XTickLabel = [];

ax1 = subaxis(2,1,2, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, TB_rel_s, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = 'Thermobaric relative error on $\\omega_{s}$-surface (OCCA), RMS = %.2f';
title(sprintf(txt ,nanrms(TB_rel_s(:))) , 'fontsize',10,'Interpreter','latex');
caxis([0, 0.1])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%%
[~, ~,TB_rel_uss2,~, ~, CB_rel_uss2] = rel_TB_CB(t_hels2xy, z_hels2xy,detail_hel_s2xy.sx, detail_hel_s2xy.sy, s2hel_s2xy, ZB, pB, TB, g.DXCvec, g.DYCsc, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

figure('Position', [0 0 1200 1000])
ax1 = subaxis(2,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, CB_rel_uss2, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = 'Cabeling relative error on $\\omega_{u.s +s^2}$-surface (OCCA), RMS = %.2f';
title(sprintf(txt ,nanrms(CB_rel_uss2(:))) , 'fontsize',10,'Interpreter','latex');
caxis([0, 0.3])
ax1.FontSize = 15;
ax1.XTickLabel = [];

ax1 = subaxis(2,1,2, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, TB_rel_uss2, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = 'Thermobaric relative error on $\\omega_{u.s+s^2}$-surface (OCCA), RMS = %.2f';
title(sprintf(txt ,nanrms(TB_rel_uss2(:))) , 'fontsize',10,'Interpreter','latex');
caxis([0, 0.3])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%%
OP.tolz = 1e-8; % tolerance of slope error
OP.dx = g.DXCvec;
OP.dy = g.DYCsc;
OP.INTERPFN = @ppc_pchip;

[grad_T_x_n_ans, grad_T_y_n_ans] = grad_n_T_fun(SB, TB, ZB, zns, OP);

grad_T_s_x = (tns_s - im1(tns_s))./g.DXCvec;
grad_T_s_y = (tns_s - jm1(tns_s))./g.DYCsc;

grad_T_x_n_ans_dely = grad_T_x_n_ans .* (g.DYGsc);
grad_T_y_n_ans_delx = grad_T_y_n_ans .* (g.DXGvec);
grad_T_n_ans_int = (grad_T_x_n_ans_dely + ip1(grad_T_x_n_ans_dely) + grad_T_y_n_ans_delx + jp1(grad_T_y_n_ans_delx)) ./ (g.RACvec);

grad_T_x_dely = grad_T_s_x .* (g.DYGsc);
grad_T_y_delx = grad_T_s_y .* (g.DXGvec);
grad_T_int = (grad_T_x_dely + ip1(grad_T_x_dely) + grad_T_y_delx + jp1(grad_T_y_delx)) ./ (g.RACvec);

ep_rel = 1000*(grad_T_n_ans_int - grad_T_int);

figure('Position', [0 0 1800 1000])
ax1 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, ep_rel, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = 'Cabeling relative error on $\\omega_{s}$-surface (OCCA), RMS = %.2d';
title(sprintf(txt ,nanrms(ep_rel(:))) , 'fontsize',10,'Interpreter','latex');
% caxis([0, 0.1])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%% relative errors on omega_s
sxTz = detailns_s.sx.* Tz_ns_s;
syTz = detailns_s.sy.* Tz_ns_s;

sxTz_dely = sxTz .* (g.DYGsc);
syTz_delx = syTz .* (g.DXGvec);
sTz_int_ns = (sxTz_dely + ip1(sxTz_dely) + syTz_delx + jp1(syTz_delx)) ./ (g.RACvec);
ep_er_ns = 1000*sTz_int_ns;
ep_rel_ns = sTz_int_ns./grad_T_int;

figure('Position', [0 0 1800 1000])
ax1 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, ep_rel_ns, land_OCCA, OPTS_FIGS)
colorbar(ax1)
txt = 'Cabeling relative error on $\\omega_{s}$-surface (OCCA), RMS = %.2f';
title(sprintf(txt ,nanrms(ep_rel_ns(:))) , 'fontsize',10,'Interpreter','latex');
caxis([0, 0.1])
ax1.FontSize = 15;
ax1.XTickLabel = [];

%% relative errors on omega_us

grad_T_s_x = (t_hels2xy - im1(t_hels2xy))./g.DXCvec;
grad_T_s_y = (t_hels2xy - jm1(t_hels2xy))./g.DYCsc;

grad_T_x_dely = grad_T_s_x .* (g.DYGsc);
grad_T_y_delx = grad_T_s_y .* (g.DXGvec);
grad_T_int = (grad_T_x_dely + ip1(grad_T_x_dely) + grad_T_y_delx + jp1(grad_T_y_delx)) ./ (g.RACvec);
Kgrad_T_int = 1000*grad_T_int;
sxTz = detailns_s2xy.sx.* Tz_hels2xy;
syTz = detailns_s2xy.sy.* Tz_hels2xy;

sxTz_dely = sxTz .* (g.DYGsc);
syTz_delx = syTz .* (g.DXGvec);
sTz_int_ns = (sxTz_dely + ip1(sxTz_dely) + syTz_delx + jp1(syTz_delx)) ./ (g.RACvec);
ep_er_ns = 1000*sTz_int_ns;
ep_rel_ns = sTz_int_ns./grad_T_int;


[TB_er,grad_p_grad_T, TB_rel,grad_T_grad_T,CB_er, CB_rel]= rel_TB_CB(t_hels2xy, z_hels2xy,detail_hel_s2xy.sx, detail_hel_s2xy.sy, s2hel_s2xy, ZB, pB, TB, g.DXCvec, g.DYCsc, g.DXCvec .* g.DYGsc, g.DXGvec .* g.DYCsc, g.RACvec);

%%
OPTS_AXES = {'Margin', .04, 'Spacing', .03};
figure('Position', [0 0 1800 1000])

ax5 = subaxis(3,2,1, OPTS_AXES{:});
hf = fig_map(ax5, g.XCvec, g.YCvec, ep_er_ns, land_OCCA, OPTS_FIGS);
colorbar(ax5)
txt = '(a) $(h^{-1} \\nabla_n \\cdot (h K \\nabla_n \\Theta) - h^{-1} \\nabla_a \\cdot (h K \\nabla_a \\Theta))$ on the $\\omega_{u.s+s^2}$-surface (OCCA), RMS = %.2d K/s';
title(sprintf(txt ,nanrms(ep_er_ns(:))) , 'fontsize',10,'Interpreter','latex');
pcol_logsigned(ax5, g.XCvec, g.YCvec, ep_er_ns, -14, -2);
ax5.FontSize = 15;
ax5.XTickLabel = [];


ax6 = subaxis(3,2,2, OPTS_AXES{:});
hf = fig_map(ax6, g.XCvec, g.YCvec, Kgrad_T_int, land_OCCA, OPTS_FIGS);
colorbar(ax6)
txt = '(b) $(h^{-1} \\nabla_a \\cdot (h K \\nabla_a \\Theta))$ on the $\\omega_{u.s+s^2}$-surface (OCCA), RMS = %.2d K/s';
title(sprintf(txt ,nanrms(Kgrad_T_int(:))) , 'fontsize',10,'Interpreter','latex');
pcol_logsigned(ax6, g.XCvec, g.YCvec, Kgrad_T_int, -14, -2);
ax6.FontSize = 15;
ax6.YTickLabel = [];
ax6.XTickLabel = [];

ax3 = subaxis(3,2,3, OPTS_AXES{:});
hf = fig_map(ax3, g.XCvec, g.YCvec, TB_er, land_OCCA, OPTS_FIGS);
colorbar(ax3)
txt = '(c) $(\\nabla_n p \\cdot \\nabla_n \\Theta - \\nabla_a p \\cdot \\nabla_a \\Theta)$ on the $\\omega_{u.s+s^2}$-surface (OCCA), RMS = %.2d $Pa K/m^2$';
title(sprintf(txt ,nanrms(TB_er(:))) , 'fontsize',10,'Interpreter','latex');
pcol_logsigned(ax3, g.XCvec, g.YCvec, TB_er, -14, -2);
ax3.FontSize = 15;
ax3.XTickLabel = [];

ax4 = subaxis(3,2,4, OPTS_AXES{:});
hf = fig_map(ax4, g.XCvec, g.YCvec, grad_p_grad_T, land_OCCA, OPTS_FIGS);
colorbar(ax4)
txt = '(d) $(\\nabla_a p \\cdot \\nabla_a \\Theta)$ on the $\\omega_{u.s+s^2}$-surface (OCCA), RMS = %.2d $Pa K/m^2$';
title(sprintf(txt ,nanrms(grad_p_grad_T(:))) , 'fontsize',10,'Interpreter','latex');
pcol_logsigned(ax4, g.XCvec, g.YCvec, grad_p_grad_T, -14, -2);
ax4.FontSize = 15;
ax4.YTickLabel = [];
ax4.XTickLabel = [];

ax1 = subaxis(3,2,5, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, CB_er, land_OCCA, OPTS_FIGS);
colorbar(ax1)
txt = '(e) $(\\nabla_n \\Theta \\cdot \\nabla_n \\Theta - \\nabla_a \\Theta \\cdot \\nabla_a \\Theta$) on the $\\omega_{u.s+s^2}$-surface (OCCA), RMS = %.2d $K^2/m^2$';
title(sprintf(txt ,nanrms(CB_er(:))) , 'fontsize',10,'Interpreter','latex');
pcol_logsigned(ax1, g.XCvec, g.YCvec, CB_er, -14, -2);
ax1.FontSize = 15;

ax2 = subaxis(3,2,6, OPTS_AXES{:});
hf = fig_map(ax2, g.XCvec, g.YCvec, grad_T_grad_T, land_OCCA, OPTS_FIGS);
colorbar(ax2)
txt = '(f) $(\\nabla_a \\Theta \\cdot \\nabla_a \\Theta)$ on the $\\omega_{u.s+s^2}$-surface (OCCA), RMS = %.2d $K^2/m^2$';
title(sprintf(txt ,nanrms(grad_T_grad_T(:))) , 'fontsize',10,'Interpreter','latex');
pcol_logsigned(ax2, g.XCvec, g.YCvec, grad_T_grad_T, -14, -2);
ax2.FontSize = 15;
ax2.YTickLabel = [];




