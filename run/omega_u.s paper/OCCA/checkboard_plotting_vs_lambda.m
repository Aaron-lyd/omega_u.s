
im1 = @(F) circshift(F, [+1 0]);
jm1 = @(F) circshift(F, [0 +1]);
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
%% omega_1
OPTS = [];
OPTS.POISSON = 0;
OPTS.INTERPFN = interpfn;

[zns, sns, tns, dns] = omega_surface(SB, TB, ZB, z_sigma, I0,[1;1], OPTS); %in x-y direction

%% Grid information
OPTS.nx = g.nx; 
OPTS.ny = g.ny;
OPTS.DXCvec = g.DXCvec;
OPTS.DXGvec = g.DXGvec;
OPTS.DYCsc = g.DYCsc;
OPTS.DYGsc = g.DYGsc;
OPTS.RACvec = g.RACvec;
%OPTS.RACvec = 1;
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
fprintf('omega_hel');
OPTS.data_cube=0;
OPTS.SHEAR = 1;
OPTS.STRAT = 1;
OPTS.PCG = 0;
OPTS.LM = 0;
OPTS.H_SIM = 0;
OPTS.CHOLESKY = 0;
OPTS.LSQLIN = 0;
OPTS.ITER = 50;
OPTS.MODE = 1;

OPTS.TK = 8e-9;
[z_hel9,s_hel9,t_hel9, ~,d_hel9] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);
OPTS.TK = 8e-5;
[z_hel5,s_hel5,t_hel5, ~,d_hel5] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);
OPTS.TK = 8e-7;
[z_hel7,s_hel7,t_hel7, ~,d_hel7] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);
OPTS.TK = 8e-11;
[z_hel11,s_hel11,t_hel11, ~,d_hel11] = omega_hel_surface(SB, TB, ZB, U, V, zns, OPTS);

%% save data
save checkboard s_hel9 t_hel9 z_hel9 d_hel9...
                s_hel7 t_hel7 z_hel7 d_hel7...
                s_hel5 t_hel5 z_hel5 d_hel5...
                s_hel11 t_hel11 z_hel11 d_hel11

%% quantify the checkerboard
dx = g.DXCvec; % dx
dy = g.DYCsc; % dy
dz_hel5_dx = squeeze((d_hel5.dz_2d(1,:,:) - im1(d_hel5.dz_2d(1,:,:))))./dx;
dz_hel5_dy = squeeze((d_hel5.dz_2d(1,:,:) - jm1(d_hel5.dz_2d(1,:,:))))./dy;
dz_hel5_rms  = rms_Cgrid(0, dz_hel5_dx, dz_hel5_dy, 0, g.RAWvec, g.RASvec);

dz_hel7_dx = squeeze((d_hel7.dz_2d(1,:,:) - im1(d_hel7.dz_2d(1,:,:))))./dx;
dz_hel7_dy = squeeze((d_hel7.dz_2d(1,:,:) - jm1(d_hel7.dz_2d(1,:,:))))./dy;
dz_hel7_rms  = rms_Cgrid(0, dz_hel7_dx, dz_hel7_dy, 0, g.RAWvec, g.RASvec);

dz_hel9_dx = squeeze((d_hel9.dz_2d(1,:,:) - im1(d_hel9.dz_2d(1,:,:))))./dx;
dz_hel9_dy = squeeze((d_hel9.dz_2d(1,:,:) - jm1(d_hel9.dz_2d(1,:,:))))./dy;
dz_hel9_rms  = rms_Cgrid(0, dz_hel9_dx, dz_hel9_dy, 0, g.RAWvec, g.RASvec);

dz_hel11_dx = squeeze((d_hel11.dz_2d(1,:,:) - im1(d_hel11.dz_2d(1,:,:))))./dx;
dz_hel11_dy = squeeze((d_hel11.dz_2d(1,:,:) - jm1(d_hel11.dz_2d(1,:,:))))./dy;
dz_hel11_rms  = rms_Cgrid(0, dz_hel11_dx, dz_hel11_dy, 0, g.RAWvec, g.RASvec);

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

[ss_hel7, df_hel7, ss_hel7_rms, df_hel7_rms, ehel7_hel, ehel_hel7_rms, dhel7] = s_ehel_df(z_hel7, SB, TB, ZB, U, V, [1;1], OPT);
[ss_hel5, df_hel5, ss_hel5_rms, df_hel5_rms, ehel5_hel, ehel_hel5_rms, dhel5] = s_ehel_df(z_hel5, SB, TB, ZB, U, V, [1;1], OPT);
[ss_hel9, df_hel9, ss_hel9_rms, df_hel9_rms, ehel9_hel, ehel_hel9_rms, dhel9] = s_ehel_df(z_hel9, SB, TB, ZB, U, V, [1;1], OPT);
[ss_hel11, df_hel11, ss_hel11_rms, df_hel11_rms, ehel11_hel, ehel_hel11_rms, dhel11] = s_ehel_df(z_hel11, SB, TB, ZB, U, V, [1;1], OPT);

%% Figure 6: Streamline
OPTS_AXES = {'Margin', .05, 'Spacing', .03};
figure('Position', [0 0 1800 1000])

x_quiver = repmat(g.XCvec, [1 160]);
y_quiver = repmat(g.YCvec, [360 1]);

OPTS_FIGS.XSH = 0;
ax1 = subaxis(2,2,1, OPTS_AXES{:});
u_hel = dhel5.u;
v_hel = dhel5.v;
hf = fig_map(ax1, g.XCvec, g.YCvec, ehel5_hel, land_OCCA, OPTS_FIGS);
pcol_logsigned(ax1, g.XCvec, g.YCvec, ehel5_hel, -11, -7);
% hold on
% hh = quivers(x_quiver, y_quiver, u_hel, v_hel,5,1,'m/s','k');
% hold off
tk5 = 8e-5;
txt = '(a) $\\lambda^2 =  %.2d$, RMS of $\\bf{u} \\cdot \\bf{s}$  = %.2d m/s, RMS of $\\nabla z_1$ = %.2d ';
title(sprintf(txt,tk5, ehel_hel5_rms, dz_hel5_rms) , 'fontsize',10,'Interpreter','latex');
ax1.FontSize = 18;
ax1.XTickLabel = [];

ax2 = subaxis(2,2,2, OPTS_AXES{:});
u_hel = dhel7.u;
v_hel = dhel7.v;
hf = fig_map(ax2, g.XCvec, g.YCvec, ehel7_hel, land_OCCA, OPTS_FIGS);
pcol_logsigned(ax2, g.XCvec, g.YCvec, ehel7_hel, -11, -7);
% hold on
% hh = quivers(x_quiver, y_quiver, u_hel, v_hel,5,1,'m/s','k');
% hold off
txt = '(b) $\\lambda^2 =  %.2d$, RMS of $\\bf{u} \\cdot \\bf{s}$  = %.2d m/s, RMS of $\\nabla z_1$ = %.2d ';
title(sprintf(txt,8e-7, ehel_hel7_rms, dz_hel7_rms) , 'fontsize',10,'Interpreter','latex');
ax2.FontSize = 18;
ax2.YTickLabel = [];
ax2.XTickLabel = [];

ax3 = subaxis(2,2,3, OPTS_AXES{:});
u_hel = dhel9.u;
v_hel = dhel9.v;
hf = fig_map(ax3, g.XCvec, g.YCvec, ehel9_hel, land_OCCA, OPTS_FIGS);
pcol_logsigned(ax3, g.XCvec, g.YCvec, ehel9_hel, -11, -7);
% hold on
% hh = quivers(x_quiver, y_quiver, u_hel, v_hel,5,1,'m/s','k');
% hold off
txt = '(c) $\\lambda^2 = %.2d$, RMS of $\\bf{u} \\cdot \\bf{s}$  = %.2d m/s, RMS of $\\nabla z_1$ = %.2d ';
title(sprintf(txt,8e-09, ehel_hel9_rms, dz_hel9_rms) , 'fontsize',10,'Interpreter','latex');
ax3.FontSize = 18;

ax4 = subaxis(2,2,4, OPTS_AXES{:});
u_hel = dhel11.u;
v_hel = dhel11.v;
hf = fig_map(ax4, g.XCvec, g.YCvec, ehel11_hel, land_OCCA, OPTS_FIGS);
pcol_logsigned(ax4, g.XCvec, g.YCvec, ehel11_hel, -11, -7);
% hold on
% hh = quivers(x_quiver, y_quiver, u_hel, v_hel,5,1,'m/s','k');
% hold off
txt = '(d) $\\lambda^2 =  %.2d$, RMS of $\\bf{u} \\cdot \\bf{s}$  = %.2d m/s, RMS of $\\nabla z_1$ = %.2d ';
title(sprintf(txt, 8e-11, ehel_hel11_rms, dz_hel11_rms) , 'fontsize',10,'Interpreter','latex');
ax4.FontSize = 18;
ax4.YTickLabel = [];


