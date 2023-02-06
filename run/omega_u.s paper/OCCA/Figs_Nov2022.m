%% s, e^hel and flux  

flat = @(x) x(:);
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
OPT.data_cube = 0;

[s_pt, df_pt, spt_rms, df_pt_rms, ehelpt, ehelpt_rms, dpt] = s_ehel_df(z_sigma, SB, TB, ZB, U, V, [1;1], OPT);
[s_ns, df_ns, sns_rms, df_ns_rms, ehelns, ehelns_rms, dns] = s_ehel_df(zns, SB, TB, ZB, U, V, [1;1], OPT);
[s_ns_s, df_ns_s, sns_s_rms, df_ns_s_rms, ehelns_s, ehelns_s_rms, dns_s] = s_ehel_df(zns_s, SB, TB, ZB, U, V, [1;1], OPT);
[ss_hel, df_hel, ss_hel_rms, df_hel_rms, ehel_hel, ehel_hel_rms, dhel] = s_ehel_df(z_hel, SB, TB, ZB, U, V, [1;1], OPT);
[ss_helTz, df_helTz, ss_helTz_rms, df_helTz_rms, ehel_helTz, ehel_helTz_rms, dhelTz] = s_ehel_df(z_helTz, SB, TB, ZB, U, V, [1;1], OPT);
[ss_helSz, df_helSz, ss_helSz_rms, df_helSz_rms, ehel_helSz, ehel_helSz_rms, dhelSz] = s_ehel_df(z_helSz, SB, TB, ZB, U, V, [1;1], OPT);
[s_ns_s2xy, df_ns_s2xy, sns_s2xy_rms, df_ns_s2xy_rms, ehelns_s2xy, ehelns_s2xy_rms, dns_s2xy] = s_ehel_df(zns_s2xy, SB, TB, ZB, U, V, [1;1], OPT);
[ss_hel_s2xy, df_hel_s2xy, ss_hel_s2xy_rms, df_hel_s2xy_rms, ehel_hel_s2xy, ehel_hel_s2xy_rms, dhel_s2xy] = s_ehel_df(z_hels2xy, SB, TB, ZB, U, V, [1;1], OPT);
%[sortho, df_ortho, sortho_rms, df_ortho_rms, ehel_ortho, ehel_ortho_rms, dortho] = s_ehel_df(z_ortho, SB, TB, ZB, U, V, [1;1], OPT);                     
%[stopo, df_topo, stopo_rms, df_topo_rms, ehel_topo, ehel_topo_rms, dtopo] = s_ehel_df(z_topo, SB, TB, ZB, U, V, [1;1], OPT);                     
[s_scv, df_scv, sscv_rms, df_scv_rms, ehel_scv, ehel_scv_rms, dscv] = s_ehel_df(zscv, SB, TB, ZB, U, V, [1;1], OPT);                     
[s_gt, df_gt, sgt_rms, df_gt_rms, ehel_gt, ehel_gt_rms, dgt] = s_ehel_df(z_gammat, SB, TB, ZB, U, V, [1;1], OPT);                     


surfacess  = {'sigma', 'omega+', 'omega_s', 'omega_hel', 'omega_helTz', 'omega_helSz', 'omega_s2xy',  'omega_hel_s2xy', 'gammaSCV', 'gammaT'}';
ehel_south_rms = [dpt.e_hel_st_rms; dns.e_hel_st_rms; dns_s.e_hel_st_rms; dhel.e_hel_st_rms; dhelTz.e_hel_st_rms; dhelSz.e_hel_st_rms; dns_s2xy.e_hel_st_rms; dhel_s2xy.e_hel_st_rms; dscv.e_hel_st_rms; dgt.e_hel_st_rms];
Df_south_rms = 1000*[dpt.s2_st_rms; dns.s2_st_rms; dns_s.s2_st_rms; dhel.s2_st_rms; dhelTz.s2_st_rms; dhelSz.s2_st_rms; dns_s2xy.s2_st_rms; dhel_s2xy.s2_st_rms; dscv.s2_st_rms; dgt.s2_st_rms];

TT = table(surfacess,ehel_south_rms,Df_south_rms)  

%% Figure 1:e^hel maps
OPTS_FIGS.LANDCOL = [1 1 1]*0; % Black
OPTS_FIGS.NANCOL = [1 1 1]* .75; % Grey
OPTS_FIGS.LATLIM = [-80, 72]; % The z(180,0) = 2000 surfaces exclude the Arctic.
land_OCCA = squeeze(isnan(SB(1,:,:)));

OPTS_plot.nx = 5;
OPTS_plot.ny = 2;
OPTS_plot.log = 1;
OPTS_plot.X = g.XCvec;
OPTS_plot.Y = g.YCvec;
OPTS_plot.land_mask = land_OCCA;
OPTS_plot.bwr = ones(1,10);
OPTS_plot.font_size = 15;
% OPTS_plot.font_size = 12;
OPTS_plot.AXES = {'Margin', .04, 'Spacing', .03};

figure('Position', [0 0 1100 1000])
OPTS_plot.n_lim = -1e-7*ones(1,10);
OPTS_plot.p_lim = 1e-7*ones(1,10);

OPTS_plot.upb = -6;
OPTS_plot.lowb = -10;

data = cat(3, ehelpt, ehelns, ehelns_s, ehel_hel, ehel_helTz, ehel_helSz, ehelns_s2xy, ehel_hel_s2xy, ehel_scv, ehel_gt);
data_rms = [ehelpt_rms, ehelns_rms, ehelns_s_rms, ehel_hel_rms, ehel_helTz_rms,...
                                                       ehel_helSz_rms, ehelns_s2xy_rms, ehel_hel_s2xy_rms, ehel_scv_rms, ehel_gt_rms];
                                                  
ehelpt_rms = rms_Cgrid(ehelpt, 0, 0, g.RACvec, 0, 0);
                                                   
quantity_txt = '$\\bf{u} \\cdot \\bf{s}$ on the ';        
rms_txt = 'RMS of $\\bf{u} \\cdot \\bf{s}$ = %.2d m/s';        
txt1 = ['(a) ', quantity_txt,'$\\sigma_{1.2}$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{\\mathbf{u} \\cdot \\mathbf{s}}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{\\mathbf{u} \\cdot \\mathbf{s}\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{\\mathbf{u} \\cdot \\mathbf{s}S_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{\\mathbf{s}^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{\\mathbf{u} \\cdot \\mathbf{s}+\\mathbf{s}^2}$-surface, ', rms_txt];
txt9 = ['(i) ', quantity_txt,'$\\gamma^{SCV}=27.98$-surface, ', rms_txt];
txt10 = ['(j) ', quantity_txt,'$\\gamma^T=27.98$-surface, ', rms_txt];

title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8,txt9,txt10);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% e_hel Tz map
ehelTz_pt = dpt.ehelTz;
ehelTz_ns = dns.ehelTz;
ehelTz_ns_s = dns_s.ehelTz;
ehelTz_hel = dhel.ehelTz;
ehelTz_helTz = dhelTz.ehelTz;
ehelTz_helSz = dhelSz.ehelTz;
ehelTz_s2xy = dns_s2xy.ehelTz;
ehelTz_hels2xy = dhel_s2xy.ehelTz;
ehelTz_scv = dscv.ehelTz;
ehelTz_gt = dgt.ehelTz;


figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -1e-9*ones(1,10);
OPTS_plot.p_lim = 1e-9*ones(1,10);
OPTS_plot.bwr = ones(1,10);

data = cat(3, ehelTz_pt, ehelTz_ns, ehelTz_ns_s, ehelTz_hel,...
                                                 ehelTz_helTz, ehelTz_helSz, ehelTz_s2xy, ehelTz_hels2xy, ehelTz_scv, ehelTz_gt);
data_rms = [dpt.ehelTz_rms,dns.ehelTz_rms,dns_s.ehelTz_rms,dhel.ehelTz_rms,dhelTz.ehelTz_rms,dhelSz.ehelTz_rms,dns_s2xy.ehelTz_rms,dhel_s2xy.ehelTz_rms,dscv.ehelTz_rms,dgt.ehelTz_rms];

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
txt9 = ['(i) ', quantity_txt,'$\\gamma^{SCV}=27.78$-surface, ', rms_txt];
txt10 = ['(j) ', quantity_txt,'$\\gamma^T=27.78$-surface, ', rms_txt];

title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8,txt9,txt10);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% e_hel Sz map
ehelSz_pt = dpt.ehelSz;
ehelSz_ns = dns.ehelSz;
ehelSz_ns_s = dns_s.ehelSz;
ehelSz_hel = dhel.ehelSz;
ehelSz_helTz = dhelTz.ehelSz;
ehelSz_helSz = dhelSz.ehelSz;
ehelSz_s2xy = dns_s2xy.ehelSz;
ehelSz_hels2xy = dhel_s2xy.ehelSz;
ehelSz_scv = dscv.ehelSz;
ehelSz_gt = dgt.ehelSz;

figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -1e-10*ones(1,10);
OPTS_plot.p_lim = 1e-10*ones(1,10);
OPTS_plot.bwr = ones(1,10);

data = cat(3, ehelSz_pt, ehelSz_ns, ehelSz_ns_s, ehelSz_hel,...
             ehelSz_helTz, ehelSz_helSz, ehelSz_s2xy, ehelSz_hels2xy, ehelSz_scv, ehelSz_gt);
data_rms = [dpt.ehelSz_rms,dns.ehelSz_rms,dns_s.ehelSz_rms,dhel.ehelSz_rms,dhelTz.ehelSz_rms,dhelSz.ehelSz_rms,dns_s2xy.ehelSz_rms,dhel_s2xy.ehelSz_rms,dscv.ehelSz_rms,dgt.ehelSz_rms];

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
txt9 = ['(i) ', quantity_txt,'$\\gamma^{SCV}=27.78$-surface, ', rms_txt];
txt10 = ['(j) ', quantity_txt,'$\\gamma^T=27.78$-surface, ', rms_txt];

title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8,txt9,txt10);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 2: Slope Error map
figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -8*ones(1,10);
OPTS_plot.p_lim = -3*ones(1,10);
OPTS_plot.bwr = zeros(1,10);

data = cat(3, log10(s_pt), log10(s_ns), log10(s_ns_s), log10(ss_hel),...
                                   log10(ss_helTz), log10(ss_helSz), log10(s_ns_s2xy), log10(ss_hel_s2xy), log10(s_scv), log10(s_gt));
data_rms = [spt_rms, sns_rms, sns_s_rms, ss_hel_rms, ss_helTz_rms,...
                                                       ss_helSz_rms, sns_s2xy_rms, ss_hel_s2xy_rms, sscv_rms, sgt_rms];
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
txt9 = ['(i) ', quantity_txt,'$\\gamma^{SCV}=27.78$-surface, ', rms_txt];
txt10 = ['(j) ', quantity_txt,'$\\gamma^T=27.78$-surface, ', rms_txt];

title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8,txt9,txt10);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 3: K s^2 map
figure('Position', [0 0 1100 1000])
OPTS_plot.n_lim = -16*ones(1,10);
OPTS_plot.p_lim = -5*ones(1,10);
OPTS_plot.bwr = zeros(1,10);

data = cat(3, log10(df_pt), log10(df_ns), log10(df_ns_s), log10(df_hel),...
                                   log10(df_helTz), log10(df_helSz), log10(df_ns_s2xy), log10(df_hel_s2xy), log10(df_scv), log10(df_gt));
data_rms = [df_pt_rms, df_ns_rms, df_ns_s_rms, df_hel_rms, df_helTz_rms,...
                                              df_helSz_rms, df_ns_s2xy_rms, df_hel_s2xy_rms, df_scv_rms, df_gt_rms];
quantity_txt = '$log_{10} (D^f)$ on the ';        
rms_txt = 'RMS = %.2d m$^2$ s$^{-1}$';        
txt1 = ['(a) ', quantity_txt,'$\\sigma_{1.2}$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{\\bf{u} \\cdot \\bf{s}S_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{\\bf{s}^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}+\\bf{s}^2}$-surface, ', rms_txt];
txt9 = ['(i) ', quantity_txt,'$\\gamma^{SCV}=27.98$-surface, ', rms_txt];
txt10 = ['(j) ', quantity_txt,'$\\gamma^T=27.98$-surface, ', rms_txt];

title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8,txt9,txt10);

fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 4: The plot of cos of omega+, omega_hel and sigma
sintheta_ns = ehelns./(abs(s_ns).*abs(dns.u_t));
sintheta_hel = ehel_hel./(abs(ss_hel).*abs(dhel.u_t));
sintheta_pt = ehelpt./(abs(s_pt).*abs(dpt.u_t));


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
OPTS_plot.n_lim = [-70,   -1800, -10, -20, -20, -20, -10, -30, -100, -100];
OPTS_plot.p_lim = [70,     0,  10,  20,  20,  20,  10,  30, 100, 100];
OPTS_plot.bwr = [1, 0, 1, 1, 1, 1, 1, 1, 1, 1];
OPTS_plot.log = 0;

data = -cat(3, (z_sigma-zns), zns,  (zns_s-zns), (z_hel-zns),...
 (z_helTz-z_hel), (z_helSz-z_hel), (zns_s2xy-zns), (z_hels2xy-z_hel), (z_sigma-zscv),(z_sigma-z_gammat));
data_rms = [nanrms(flat(z_sigma-zns)), rms(-zns(~isnan(zns))), ...
            nanrms(flat(zns_s-zns)),    nanrms(flat(z_hel-zns)),...
            nanrms(flat(z_helTz-z_hel)),  nanrms(flat(z_helSz-z_hel)),...
            nanrms(flat(zns_s2xy-zns)), nanrms(flat(z_hels2xy-z_hel)),...
            nanrms(flat(zscv-zns)), nanrms(flat(z_gammat-zns))];
       
txt1 = '(a) z[$\\sigma_{0.75}]$ - z[$\\omega_+$], RMS = %6.2f m ';
txt2 = '(b) z[$\\omega_+$], RMS = %.2d m';
txt3 = '(c) z[$\\omega_s]$ - z[$\\omega_+$], RMS = %6.2f m';
txt4 = '(d) z[$\\omega_{\\bf{u} \\cdot \\bf{s}}]$ - z[$\\omega_+$], RMS = %6.2f m';
txt5 = '(e) z[$\\omega_{\\bf{u} \\cdot \\bf{s}\\Theta_z}]$ - z[$\\omega_{\\bf{u} \\cdot \\bf{s}}$], RMS = %6.2f m';
txt6 = '(f) z[$\\omega_{\\bf{u} \\cdot \\bf{s}S_z}]$ - z[$\\omega_{\\bf{u} \\cdot \\bf{s}}$], RMS = %6.2f m';
txt7 = '(g) z[$\\omega_{\\bf{s}^2}]$ - z[$\\omega_+$], RMS = %6.2f m';
txt8 = '(h) z[$\\omega_{\\bf{u} \\cdot \\bf{s}+s^2}]$ - z[$\\omega_{\\bf{u} \\cdot \\bf{s}}$], RMS = %6.2f m';
txt9 = '(i) z[$\\gamma^{SCV}$] - z[$\\omega_+$], RMS = %6.2f m';
txt10 = '(j) z[$\\gamma^{T}$] - z[$\\omega_+$], RMS = %6.2f m';

title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8,txt9,txt10);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 6: Streamline
OPTS_AXES = {'Margin', .08, 'Spacing', .08};
figure('Position', [0 0 1800 1000])

x_quiver = repmat(g.XCvec, [1 160]);
y_quiver = repmat(g.YCvec, [360 1]);
u_hel = dhel.u;
v_hel = dhel.v;
shift = [20 0];

OPTS_FIGS.XSH = 0;
ax1 = subaxis(1,1,1, OPTS_AXES{:});
hf = fig_map(ax1, g.XCvec, g.YCvec, ehel_hel, land_OCCA, OPTS_FIGS);
% colorbar(ax1)
% caxis([-5e-9, 5e-9])
% colormap(ax1, bluewhitered), colorbar
pcol_logsigned(ax1, g.XCvec, g.YCvec, ehel_hel, -10, -7);
hold on
% quiver(x_quiver, y_quiver, u_hel, v_hel,5, 'k')
hh = quivers(x_quiver, y_quiver, u_hel, v_hel,5,1,'m/s','k');
hold off

txt = '$\\bf{u} \\cdot \\bf{s}$ on the $\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface with streamlines';
title(sprintf(txt) , 'fontsize',10,'Interpreter','latex');
ax1.FontSize = 20;

%% Figure 9: The sum of the advection and diffusion
ad_df_pt = dpt.ad_df;
ad_df_ns = dns.ad_df;
ad_df_ns_s = dns_s.ad_df;
ad_df_hel = dhel.ad_df;
ad_df_helTz = dhelTz.ad_df;
ad_df_helSz = dhelSz.ad_df;
ad_df_s2xy = dns_s2xy.ad_df;
ad_df_hel_s2xy = dhel_s2xy.ad_df;
ad_df_scv = dscv.ad_df;
ad_df_gt = dgt.ad_df;

figure('Position', [0 0 1200 1000])
OPTS_plot.n_lim = -16*ones(1,10);
OPTS_plot.p_lim = -5*ones(1,10);
OPTS_plot.bwr = zeros(1,10);
OPTS_plot.font_size = 13;
data = cat(3, ad_df_pt, ad_df_ns, ad_df_ns_s, ad_df_hel,...
                                          ad_df_helTz, ad_df_helSz, ad_df_s2xy, ad_df_hel_s2xy, ad_df_scv, ad_df_gt);
data = log10(data);
        
data_rms = [dpt.ad_df_rms, dns.ad_df_rms, dns_s.ad_df_rms, dhel.ad_df_rms, dhelTz.ad_df_rms,...
                                       dhelSz.ad_df_rms,dns_s2xy.ad_df_rms,dhel_s2xy.ad_df_rms, dscv.ad_df_rms, dgt.ad_df_rms];
        
        
quantity_txt = '$log_{10}(|\\bf{u} \\cdot \\bf{s}|$ + $10^{-3} (m^{-1})*D^f)$ on the ';        
rms_txt = 'RMS  = %.2d m/s';        
txt1 = ['(a) ', quantity_txt,'$\\sigma_{0.75}$-surface, ', rms_txt];
txt2 = ['(b) ', quantity_txt,'$\\omega_+$-surface, ', rms_txt];
txt3 = ['(c) ', quantity_txt,'$\\omega_s$-surface, ', rms_txt];
txt4 = ['(d) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}}$-surface, ', rms_txt];
txt5 = ['(e) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}\\Theta_z}$-surface, ', rms_txt];
txt6 = ['(f) ', quantity_txt,'$$\\omega_{\\bf{u} \\cdot \\bf{s}S_z}$-surface, ', rms_txt];
txt7 = ['(g) ', quantity_txt,'$\\omega_{\\bf{s}^2}$-surface, ', rms_txt];
txt8 = ['(h) ', quantity_txt,'$\\omega_{\\bf{u} \\cdot \\bf{s}+ \\bf{s}^2}$-surface, ', rms_txt];
txt9 = ['(i) ', quantity_txt,'$\\gamma^{SCV}=27.78$-surface, ', rms_txt];
txt10 = ['(j) ', quantity_txt,'$\\gamma^T=27.78$-surface, ', rms_txt];

title_text = char(txt1,txt2,txt3,txt4,txt5,txt6,txt7,txt8,txt9,txt10);
fig_hf = fig_map_plotting(data, data_rms, title_text, OPTS_plot, OPTS_FIGS);

%% Figure 10: The scatter of e^hel and s^2 for omega_hel and omega_hel+s^2
OPTS_AXES = {'Margin', .06, 'Spacing', .09};
figure('Position', [0 0 1500 1000])

ax1 = subaxis(3,2,1, OPTS_AXES{:});
good = ~isnan(ehelpt) & ~isnan(df_pt);
xvec = ehelpt(good);
yvec = df_pt(good);
hist2D(ax1,log10(abs(xvec)), log10(abs(yvec*1e-3)),'nbins',400, 'RescaleFcn', @(x) sqrt(x))
plot(log10(ehelpt_rms), log10(dpt.s2_rms), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
patch_x = [-12 log10(3e-9) log10(3e-9) -12];
patch_y = [log10(1e-8) log10(1e-8) -20 -20];
blue = [0.5843 0.8157 0.9882];
patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
grid on
colorbar
xlim([-12, -4])
ylim([-20, -5])
% refline(1,0)
xline(log10(3e-9),'-',{'Acceptable','Limit'});
yline(log10(1e-8),'-',{'Acceptable','Limit'});
txt_x = '$log_{10} (\\ |\\bf{u} \\cdot \\bf{s}\\ |)$ [m/s], RMS = %.2d m/s';
xlabel(sprintf(txt_x ,ehelpt_rms) , 'fontsize',10,'Interpreter','latex');
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';
ylabel(sprintf(txt_y ,dpt.s2_rms) , 'fontsize',10,'Interpreter','latex');
title('(a) $\sigma_{0.75}$-surface', 'fontsize',10,'Interpreter','latex')
set(gca,'fontsize', 16);

ax2 = subaxis(3,2,2, OPTS_AXES{:});
good = ~isnan(ehelns) & ~isnan(df_ns);
xvec = ehelns(good);
yvec = df_ns(good);
hist2D(ax2,log10(abs(xvec)), log10(abs(yvec*1e-3)),'nbins',400, 'RescaleFcn', @(x) sqrt(x))
plot(log10(ehelns_rms), log10(dns.s2_rms), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
patch_x = [-12 log10(3e-9) log10(3e-9) -12];
patch_y = [log10(1e-8) log10(1e-8) -20 -20];
patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
grid on
colorbar
xlim([-12, -4])
ylim([-20, -5])
% refline(1,0)
xline(log10(3e-9),'-',{'Acceptable','Limit'});
yline(log10(1e-8),'-',{'Acceptable','Limit'});
txt_x = '$log_{10} (\\ |\\bf{u} \\cdot \\bf{s}\\ |)$[m/s], RMS = %.2d m/s';
xlabel(sprintf(txt_x ,ehelns_rms) , 'fontsize',10,'Interpreter','latex');
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';
ylabel(sprintf(txt_y ,dns.s2_rms) , 'fontsize',10,'Interpreter','latex');
title('(b) $\omega_{+}$-surface', 'fontsize',10,'Interpreter','latex')
set(gca,'fontsize', 16);

ax3 = subaxis(3,2,3, OPTS_AXES{:});
good = ~isnan(ehel_hel) & ~isnan(df_hel);
xvec = ehel_hel(good);
yvec = df_hel(good);
hist2D(ax3,log10(abs(xvec)), log10(abs(yvec*1e-3)),'nbins',400, 'RescaleFcn', @(x) sqrt(x))
plot(log10(ehel_hel_rms), log10(dhel.s2_rms), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
blue = [ 0.5843 0.8157 0.9882];
patch_x = [-12 log10(3e-9) log10(3e-9) -12];
patch_y = [log10(1e-8) log10(1e-8) -20 -20];
patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
grid on
colorbar
xlim([-12, -4])
ylim([-20, -5])
% refline(1,0)
xline(log10(3e-9),'-',{'Acceptable','Limit'});
yline(log10(1e-8),'-',{'Acceptable','Limit'});
txt_x = '$log_{10} (\\ |\\bf{u} \\cdot \\bf{s}\\ |)$ [m/s], RMS = %.2d m/s';
xlabel(sprintf(txt_x ,ehel_hel_rms) , 'fontsize',10,'Interpreter','latex');
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';
ylabel(sprintf(txt_y ,dhel.s2_rms) , 'fontsize',10,'Interpreter','latex');
title('(c) $\omega_{\bf{u} \cdot \bf{s}}$-surface', 'fontsize',10,'Interpreter','latex')
set(gca,'fontsize', 16);

ax4 = subaxis(3,2,4, OPTS_AXES{:});
good = ~isnan(ehel_hel_s2xy) & ~isnan(df_hel_s2xy);
xvec = ehel_hel_s2xy(good);
yvec = df_hel_s2xy(good);
hist2D(ax4,log10(abs(xvec)), log10(abs(yvec*1e-3)),'nbins',400, 'RescaleFcn', @(x) sqrt(x))
plot(log10(ehel_hel_s2xy_rms), log10(dhel_s2xy.s2_rms), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
patch_x = [-12 log10(3e-9) log10(3e-9) -12];
patch_y = [log10(1e-8) log10(1e-8) -20 -20];
patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
grid on
colorbar
xlim([-12, -4])
ylim([-20, -5])
% refline(1,0)
xline(log10(3e-9),'-',{'Acceptable','Limit'});
yline(log10(1e-8),'-',{'Acceptable','Limit'});
txt_x = '$log_{10} (\\ |\\bf{u} \\cdot \\bf{s}\\ |)$ [m/s], RMS = %.2d m/s';
xlabel(sprintf(txt_x ,ehel_hel_s2xy_rms) , 'fontsize',10,'Interpreter','latex');
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';
ylabel(sprintf(txt_y ,dhel_s2xy.s2_rms) , 'fontsize',10,'Interpreter','latex');
title('(d) $\omega_{\bf{u} \cdot \bf{s}+\bf{s}^2}$-surface', 'fontsize',10,'Interpreter','latex')
set(gca,'fontsize', 16);

ax5 = subaxis(3,2,5, OPTS_AXES{:});
good = ~isnan(ehel_scv) & ~isnan(df_scv);
xvec = ehel_scv(good);
yvec = df_scv(good);
hist2D(ax5,log10(abs(xvec)), log10(abs(yvec*1e-3)),'nbins',400, 'RescaleFcn', @(x) sqrt(x))
plot(log10(ehel_scv_rms), log10(dscv.s2_rms), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
patch_x = [-12 log10(3e-9) log10(3e-9) -12];
patch_y = [log10(1e-8) log10(1e-8) -20 -20];
patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
grid on
colorbar
xlim([-12, -4])
ylim([-20, -5])
% refline(1,0)
xline(log10(3e-9),'-',{'Acceptable','Limit'});
yline(log10(1e-8),'-',{'Acceptable','Limit'});
txt_x = '$log_{10} (\\ |\\bf{u} \\cdot \\bf{s}\\ |)$ [m/s], RMS = %.2d m/s';
xlabel(sprintf(txt_x ,ehel_scv_rms) , 'fontsize',10,'Interpreter','latex');
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';
ylabel(sprintf(txt_y ,dscv.s2_rms) , 'fontsize',10,'Interpreter','latex');
title('(e) $\gamma^{SCV}=27.78$-surface', 'fontsize',10,'Interpreter','latex')
set(gca,'fontsize', 16);

ax6 = subaxis(3,2,6, OPTS_AXES{:});
good = ~isnan(ehel_gt) & ~isnan(df_gt);
xvec = ehel_gt(good);
yvec = df_gt(good);
hist2D(ax6,log10(abs(xvec)), log10(abs(yvec*1e-3)),'nbins',400, 'RescaleFcn', @(x) sqrt(x))
plot(log10(ehel_gt_rms), log10(dgt.s2_rms), 'Marker', '+', 'Color', [0 0 0], 'MarkerSize', 20, 'LineWidth', 3)
patch_x = [-12 log10(3e-9) log10(3e-9) -12];
patch_y = [log10(1e-8) log10(1e-8) -20 -20];
patch(patch_x, patch_y, -1 + zeros(length(patch_x),1), blue, 'EdgeColor', 'none')
grid on
colorbar
xlim([-12, -4])
ylim([-20, -5])
% refline(1,0)
xline(log10(3e-9),'-',{'Acceptable','Limit'});
yline(log10(1e-8),'-',{'Acceptable','Limit'});
txt_x = '$log_{10} (\\ |\\bf{u} \\cdot \\bf{s}\\ |)$ [m/s], RMS = %.2d m/s';
xlabel(sprintf(txt_x ,ehel_gt_rms) , 'fontsize',10,'Interpreter','latex');
txt_y = '$log_{10} (\\bf{s}^2)$, RMS = %.2d';
ylabel(sprintf(txt_y ,dgt.s2_rms) , 'fontsize',10,'Interpreter','latex');
title('(f) $\gamma^T=27.78$-surface', 'fontsize',10,'Interpreter','latex')
set(gca,'fontsize', 16);