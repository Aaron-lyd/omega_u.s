%% Figure 1:e^hel maps
OPTS_FIGS.LANDCOL = [1 1 1]*0; % Black
OPTS_FIGS.NANCOL = [1 1 1]* .75; % Grey
OPTS_FIGS.LATLIM = [-80, 72]; % The z(180,0) = 2000 surfaces exclude the Arctic.
land_OCCA = squeeze(isnan(SB(1,:,:)));

OPTS_plot.nx = 4;
OPTS_plot.ny = 2;
OPTS_plot.log = 1;
OPTS_plot.X = g.XCvec;
OPTS_plot.Y = g.YCvec;
OPTS_plot.land_mask = land_OCCA;
OPTS_plot.bwr = ones(1,8);
OPTS_plot.font_size = 15;
OPTS_plot.AXES = {'Margin', .04, 'Spacing', .03};

figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -1e-7*ones(1,8);
OPTS_plot.p_lim = 1e-7*ones(1,8);

OPTS_plot.upb = -6;
OPTS_plot.lowb = -10;

data = cat(3, ehelpt, ehelns, ehelns_s, ehel_hel, ehel_helTz, ehel_helSz, ehelns_s2xy, ehel_hel_s2xy);
data_rms = [ehelpt_rms, ehelns_rms, ehelns_s_rms, ehel_hel_rms, ehel_helTz_rms,...
                                                       ehel_helSz_rms, ehelns_s2xy_rms, ehel_hel_s2xy_rms];

quantity_txt = '$\\bf{u} \\cdot \\bf{s}$ on the ';        
rms_txt = 'RMS of $\\bf{u} \\cdot \\bf{s}$ = %.2d m/s';        
txt1 = ['(a) ', quantity_txt,'$\\sigma_{0.75}$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{\\mathbf{u} \\cdot \\mathbf{s}}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{\\mathbf{u} \\cdot \\mathbf{s}\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{\\mathbf{u} \\cdot \\mathbf{s}S_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{\\mathbf{s}^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{\\mathbf{u} \\cdot \\mathbf{s}+\\mathbf{s}^2}$-surface, ', rms_txt];

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
quantity_txt = '$\\mathbf{u} \\cdot \\mathbf{s}\\Theta_z$ on the ';        
rms_txt = 'RMS of $\\Theta_z e^{hel}$ = %.2d $^\\circ$C/s';        
txt1 = ['(a) ', quantity_txt,'$\\sigma_{0.75}$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{hel}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{hel\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{helS_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{s^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{hel+s^2}$-surface, ', rms_txt];
title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Tz map
figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -3e-2*ones(1,8);
OPTS_plot.p_lim = 3e-2*ones(1,8);
OPTS_plot.bwr = ones(1,8);

data = cat(3, Tz_pt, Tz_ns, Tz_ns_s, Tz_hel, Tz_helTz, Tz_helSz, Tz_s2xy, Tz_hels2xy);
data_rms = [root_mean_square(Tz_pt(:)),    root_mean_square(Tz_ns(:)),...
            root_mean_square(Tz_ns_s(:)),  root_mean_square(Tz_hel(:)),...
            root_mean_square(Tz_helTz(:)), root_mean_square(Tz_helSz(:)),...
            root_mean_square(Tz_s2xy(:)),  root_mean_square(Tz_hels2xy(:))];
quantity_txt = '$\\Theta_z e^{hel}$ on ';        
rms_txt = 'RMS of $\\Theta_z e^{hel}$ = %.2d $^\\circ$C/s';        
txt1 = ['(a) ', quantity_txt,'$\\sigma_{0.75}$-surface, ', rms_txt];
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
txt1 = ['(a) ', quantity_txt,'$\\sigma_{0.75}$-surface, ', rms_txt];
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
quantity_txt = '$log_{10} (|\\bf{s}|)$ on the ';        
rms_txt = 'RMS = %.2d';        
txt1 = ['(a) ', quantity_txt,'$\\sigma_{0.75}$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{\\bf{u} \\cdot \\bf{s} S_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{\\bf{s}^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}+\\bf{s}^2}$-surface, ', rms_txt];
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
quantity_txt = '$log_{10} (D^f)$ on the ';        
rms_txt = 'RMS = %.2d m$^2$ s$^{-1}$';        
txt1 = ['(a) ', quantity_txt,'$\\sigma_{0.75}$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{\\bf{u} \\cdot \\bf{s}S_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{\\bf{s}^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}+\\bf{s}^2}$-surface, ', rms_txt];
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
txt = '(a) $\\mathbf{u} \\cdot \\mathbf{s}/|\\mathbf{u}||\\mathbf{s}|$ on the $\\omega_+$-surface';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
caxis([-1, 1])
ax1.FontSize = 15;
ax1.XTickLabel = [];

ax2 = subaxis(3,1,2, OPTS_AXES{:});
hf = fig_map(ax2, g.XCvec, g.YCvec, sintheta_hel , land_OCCA, OPTS_FIGS);
colorbar(ax2)
txt = '(b) $\\mathbf{u} \\cdot \\mathbf{s}/|\\mathbf{u}||\\mathbf{s}|$ on the $\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface';
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
txt = '(c) $\\mathbf{u} \\cdot \\mathbf{s}/|\\mathbf{u}||\\mathbf{s}|$ on the $\\sigma_{0.75}$-surface';
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
OPTS_plot.log = 0;

data = cat(3, zns, (z_sigma-zns), (zns_s-zns), (z_hel-zns),...
                                       (z_helTz-z_hel), (z_helSz-z_hel), (zns_s2xy-zns), (z_hels2xy-z_hel));
data_rms = [root_mean_square(zns),          root_mean_square(z_sigma-zns),...
            root_mean_square(zns_s-zns),    root_mean_square(z_hel-zns),...
            root_mean_square(z_helTz-z_hel),  root_mean_square(z_helSz-z_hel),...
            root_mean_square(zns_s2xy-zns), root_mean_square(z_hels2xy-z_hel)];
       
txt1 = '(a) z[$\\omega_+$], RMS = %.2d m';
txt2 = '(b) z[$\\sigma_{0.75}]$ - z[$\\omega_+$], RMS = %6.2f m';
txt3 = '(c) z[$\\omega_s]$ - z[$\\omega_+$], RMS = %6.2f m';
txt4 = '(d) z[$\\omega_{\\bf{u} \\cdot \\bf{s}}]$ - z[$\\omega_+$], RMS = %6.2f m';
txt5 = '(e) z[$\\omega_{\\bf{u} \\cdot \\bf{s}\\Theta_z}]$ - z[$\\omega_{\\bf{u} \\cdot \\bf{s}}$], RMS = %6.2f m';
txt6 = '(f) z[$\\omega_{\\bf{u} \\cdot \\bf{s}S_z}]$ - z[$\\omega_{\\bf{u} \\cdot \\bf{s}}$], RMS = %6.2f m';
txt7 = '(g) z[$\\omega_{\\bf{s}^2}]$ - z[$\\omega_+$], RMS = %6.2f m';
txt8 = '(h) z[$\\omega_{\\bf{u} \\cdot \\bf{s}+s^2}]$ - z[$\\omega_{\\bf{u} \\cdot \\bf{s}}$], RMS = %6.2f m';
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
% colorbar(ax1)
% caxis([-5e-9, 5e-9])
% colormap(ax1, bluewhitered), colorbar
pcol_logsigned(ax1, g.XCvec, g.YCvec, ehel_hel, -10, -7);
hold on
quiver(x_quiver, y_quiver, u_hel, v_hel,5, 'k')
hold off

txt = '$\\bf{u} \\cdot \\bf{s}$ on the $\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface with streamlines';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
ax1.FontSize = 20;

%% Figure 9: The sum of the advection and diffusion
ad_df_ns       = 1e-3*df_ns + ehelns;
ad_df_ns_s     = 1e-3*df_ns_s+ ehelns_s;
ad_df_pt       = 1e-3*df_pt + ehelpt;
ad_df_hel      = 1e-3*df_hel + ehel_hel;
ad_df_helTz      = 1e-3*df_helTz + ehel_helTz;
ad_df_helSz      = 1e-3*df_helSz + ehel_helSz;
ad_df_s2xy     = 1e-3*df_s2xy + ehelns_s2xy;
ad_df_hel_s2xy = 1e-3*df_hel_s2xy + ehel_hel_s2xy;

figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -1e-7*ones(1,8);
OPTS_plot.p_lim = 1e-7*ones(1,8);
OPTS_plot.bwr = ones(1,8);
OPTS_plot.upb = -6;
OPTS_plot.lowb = -10;
OPTS_plot.log = 1;

data = cat(3, ad_df_pt, ad_df_ns, ad_df_ns_s, ad_df_hel,...
                                          ad_df_helTz, ad_df_helSz, ad_df_s2xy, ad_df_hel_s2xy);
data_rms = [root_mean_square(ad_df_pt(:)),    root_mean_square(ad_df_ns(:)),...
            root_mean_square(ad_df_ns_s(:)),  root_mean_square(ad_df_hel(:)),...
            root_mean_square(ad_df_helTz(:)), root_mean_square(ad_df_helSz(:)),...
            root_mean_square(ad_df_s2xy(:)),  root_mean_square(ad_df_hel_s2xy(:))];
quantity_txt = '$\\bf{u} \\cdot \\bf{s}$ + $10^{-3} (m^{-1})*D^f$ on the ';        
rms_txt = 'RMS  = %.2d m/s';        
txt1 = ['(a) ', quantity_txt,'$\\sigma_{0.75}$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{\\bf{u} \\cdot \\bf{s}S_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{\\bf{s}^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}+ \\bf{s}^2}$-surface, ', rms_txt];
title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 10: The scatter of e^hel and s^2 for omega_hel and omega_hel+s^2
OPTS_AXES = {'Margin', .06, 'Spacing', .09};
figure('Position', [0 0 1800 1000])

ax1 = subaxis(2,2,1, OPTS_AXES{:});
good = ~isnan(ehelpt) & ~isnan(df_pt);
xvec = ehelpt(good);
yvec = df_pt(good);
hist2D(ax1,log10(abs(xvec)), log10(abs(yvec*1e-3)),'nbins',400, 'RescaleFcn', @(x) sqrt(x))
plot(log10(ehelpt_rms), log10(s2pt_rms), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
patch_x = [-12 log10(3e-9) log10(3e-9) -12];
patch_y = [log10(1e-8) log10(1e-8) -20 -20];
blue = [0.5843 0.8157 0.9882];
patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
grid on
colorbar
xlim([-12, -4])
ylim([-20, -5])
refline(1,0)
xline(log10(3e-9),'-',{'Acceptable','Limit'});
yline(log10(1e-8),'-',{'Acceptable','Limit'});
txt_x = '$log_{10} (\\ |\\bf{u} \\cdot \\bf{s}\\ |)$ [m/s], RMS = %.2d m/s';
xlabel(sprintf(txt_x ,ehelpt_rms) , 'fontsize',10,'Interpreter','latex');
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';
ylabel(sprintf(txt_y ,s2pt_rms) , 'fontsize',10,'Interpreter','latex');
title('(a) $\sigma_{0.75}$-surface', 'fontsize',10,'Interpreter','latex')
set(gca,'fontsize', 18);

ax2 = subaxis(2,2,2, OPTS_AXES{:});
good = ~isnan(ehelns) & ~isnan(df_ns);
xvec = ehelns(good);
yvec = df_ns(good);
hist2D(ax2,log10(abs(xvec)), log10(abs(yvec*1e-3)),'nbins',400, 'RescaleFcn', @(x) sqrt(x))
plot(log10(ehelns_rms), log10(s2ns_rms), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
patch_x = [-12 log10(3e-9) log10(3e-9) -12];
patch_y = [log10(1e-8) log10(1e-8) -20 -20];
patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
grid on
colorbar
xlim([-12, -4])
ylim([-20, -5])
refline(1,0)
xline(log10(3e-9),'-',{'Acceptable','Limit'});
yline(log10(1e-8),'-',{'Acceptable','Limit'});
txt_x = '$log_{10} (\\ |\\bf{u} \\cdot \\bf{s}\\ |)$[m/s], RMS = %.2d m/s';
xlabel(sprintf(txt_x ,ehelns_rms) , 'fontsize',10,'Interpreter','latex');
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';
ylabel(sprintf(txt_y ,s2ns_rms) , 'fontsize',10,'Interpreter','latex');
title('(b) $\omega_{+}$-surface', 'fontsize',10,'Interpreter','latex')
set(gca,'fontsize', 18);

ax3 = subaxis(2,2,3, OPTS_AXES{:});
good = ~isnan(ehel_hel) & ~isnan(df_hel);
xvec = ehel_hel(good);
yvec = df_hel(good);
hist2D(ax3,log10(abs(xvec)), log10(abs(yvec*1e-3)),'nbins',400, 'RescaleFcn', @(x) sqrt(x))
plot(log10(ehel_hel_rms), log10(s2hel_rms), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
blue = [ 0.5843 0.8157 0.9882];
patch_x = [-12 log10(3e-9) log10(3e-9) -12];
patch_y = [log10(1e-8) log10(1e-8) -20 -20];
patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
grid on
colorbar
xlim([-12, -4])
ylim([-20, -5])
refline(1,0)
xline(log10(3e-9),'-',{'Acceptable','Limit'});
yline(log10(1e-8),'-',{'Acceptable','Limit'});
txt_x = '$log_{10} (\\ |\\bf{u} \\cdot \\bf{s}\\ |)$ [m/s], RMS = %.2d m/s';
xlabel(sprintf(txt_x ,ehel_hel_rms) , 'fontsize',10,'Interpreter','latex');
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';
ylabel(sprintf(txt_y ,s2hel_rms) , 'fontsize',10,'Interpreter','latex');
title('(c) $\omega_{\bf{u} \cdot \bf{s}}$-surface', 'fontsize',10,'Interpreter','latex')
set(gca,'fontsize', 18);

ax4 = subaxis(2,2,4, OPTS_AXES{:});
good = ~isnan(ehel_hel_s2xy) & ~isnan(df_hel_s2xy);
xvec = ehel_hel_s2xy(good);
yvec = df_hel_s2xy(good);
hist2D(ax4,log10(abs(xvec)), log10(abs(yvec*1e-3)),'nbins',400, 'RescaleFcn', @(x) sqrt(x))
plot(log10(ehel_hel_s2xy_rms), log10(s2hel_s2xy_rms), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
patch_x = [-12 log10(3e-9) log10(3e-9) -12];
patch_y = [log10(1e-8) log10(1e-8) -20 -20];
patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
grid on
colorbar
xlim([-12, -4])
ylim([-20, -5])
refline(1,0)
xline(log10(3e-9),'-',{'Acceptable','Limit'});
yline(log10(1e-8),'-',{'Acceptable','Limit'});
txt_x = '$log_{10} (\\ |\\bf{u} \\cdot \\bf{s}\\ |)$ [m/s], RMS = %.2d m/s';
xlabel(sprintf(txt_x ,ehel_hel_s2xy_rms) , 'fontsize',10,'Interpreter','latex');
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';
ylabel(sprintf(txt_y ,s2hel_s2xy_rms) , 'fontsize',10,'Interpreter','latex');
title('(d) $\omega_{\bf{u} \cdot \bf{s}+\bf{s}^2}$-surface', 'fontsize',10,'Interpreter','latex')
set(gca,'fontsize', 18);
