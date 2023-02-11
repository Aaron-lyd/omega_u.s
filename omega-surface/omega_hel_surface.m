function [zns, sns, tns, esp, d] = omega_hel_surface(S, T, Z, U, V, z, OPTS)
% Create a series of omega surfaces, users can select the mode (in the the options)
% to choose what surface to create, OPTS.MODE =
% 1. omega_hel surface, minimizes the L2 norm of the spurious diasurface velocity; 
% 2. omega_s surface, minimizes the L2 norm of the slope error;
% 3. omega_1.5 surface, minimizes the L2 norm of e^hel_x and e^hel_y on the u and v cell
% of the c-grid;
% 4. omega_s^2, minimizes the L2 norm of the squared slope error at the tracer cell;
% 5. combine omega_hel (MODE=1) and omega_s^2 (MODE=4);
% 6. omega_s2xy, minimizes the L2 norm of the squared slope error at the u and v cell;
% 7. combine omega_hel (MODE=1) and omega_1.5 (MODE=3);
% 8. combine omega_hel (MODE=1) and omega_s2xy (MODE=6);
% 9. omega_hel weighted by Tz
% 10.omega_hel weighted by Sz
%
% [zns, sns, tns, ehel, d] = omega_hel_surface(S, T, Z, z, ref_cast, WRAP)
% returns the pressure (or depth) zns, practical / Absolute salinity sns, and
% potential / Conservative temperature tns, spurious dia-surface velocity
% esp, diagnostic terms d
% on an omega surface, 
% initialized from an approximately neutral surface of (input) pressure (or depth) z,
% in an ocean whose practical / Absolute salinity and potential /
% Conservative temperature are S and T located at datasites where the
% pressure (or depth) is Z.  

% ... = omega_surface(..., OPTS)
% specifies algorithmic parameters (see "Options" below for details).  
%
%
% --- Input:
% S [nk, ni, nj]: practical / Absolute Salinity
% T [nk, ni, nj]: potential / Conservative Temperature
% Z [nk, ni, nj] or [nk, 1]: pressure (or depth)
% U [nk, ni, nj]: zonal velocity
% V [nk, ni, nj]: meridional velocity
% z     [ni, nj]: pressure (or depth) on initial surface
%
% OPTS [struct]: options (see "Options" below)
%
%
% --- Output:
% zns [ni, nj]: pressure (or depth) on omega surface
% sns [ni, nj]: practical / Absolute salinity on omega surface
% tns [ni, nj]: potential / Conservative temperature on omega surface
% esp [ni, nj]: e^hel on omega surface
% d   [struct]: diagnostics such as clock time and norms of slope
%                 errors.  See code for info. Programmable as needed.
%
% Note: physical units of S, T, Z, and z are determined by eos.m.
%
%
% --- Equation of State:
% The MATLAB path* must contain two functions, eos.m and eos_s_t.m. Both
% accept 3 inputs: S, T, and P. eos(S, T, P) returns the specific volume
% [m^3 kg^-1] or the in-situ density [kg m^-3]. eos_s_t(S, T, P) returns,
% as its two outputs, the partial derivatives of eos with respect to S and
% T.
% *Note: It is not sufficient to simply have these eos functions in the
% current working directory, because the compiled MEX functions will not be
% able to find them there.  They must be in the MATLAB path.  If they are
% in the current working directory, use `addpath(pwd)` to add the current
% working directory to the top of MATLAB's path.
%
% For a non-Boussinesq ocean, p and P are pressure [dbar].
%
% For a Boussinesq ocean, p and P are actually depth [m].  It is essential
% that these, like pressure, are positive and increasing down.
%
% Various equation of state functions are found in ../lib/eos/.  Simply
% copy the desired functions to another location in the MATLAB path (such
% as this directory) and rename them eos.m and eos_s_t.m.  Note, the
% Boussinesq equation of state is often (but not always) just the regular
% equation of state but using a hydrostatic pressure (10^-4 * grav * rho_c
% * z) where grav [m s^-2] is the gravitational acceleration, rho_c [kg
% m^-3] is the Boussinesq reference density, and z [m, positive] is the
% depth. In such a case, simply make new eos.m and eos_p.m functions that
% accept depth as the third input by modifying the original functions that
% take pressure; this involves hard-coding the gravitational acceleration
% and Boussinesq reference density into the function.  An example of a
% Boussinesq eos.m and eos_s_t.m are given for the densjmd95 equation of
% state, in ../lib/eos/eoscg_densjmd95_bsq.m and
% ../lib/eos/eoscg_densjmd95_bsq_s_t.m.  Finally, note that eos.m and
% eos_s_t.m must be compatible with MATLAB's code generation (codegen),
% which may entail eliminating input checks and/or expansion of input
% variables (MATLAB's automatic expansion now handles this).
%
%
% --- Options:
% OPTS is a struct containing the following fields. 
%   ref_cast [1, 1] or [2, 1] : linear index or 2D index to the reference cast
%   WRAP [2 element array]: determines which dimensions are treated periodic
%                         [logical].  Set WRAP(i) to true when periodic in 
%                         the i'th lateral dimension(i=1,2).
%   FILE_ID [1, 1]: 1 to write any output to MATLAB terminal, or a file
%       identifier as returned by fopen() to write to a file. Default: 1.
%   FIGS_SHOW [scalar]: true to show figures of specific volume adjustment
%       during computation. Default: false.
%   ITER_MAX [1, 1]: maximum number of iterations. Default: 10
%   INTERPFN [function handle]: vertical interpolation function, used to
%       evaluate Sppc and Tppc if those are not provided.  Default:
%       INTERPFN = @ppc_linterp.
%   ML []: do not remove the mixed layer (default)
%   ML [struct]: calculate the mixed layer using these parameters in mixed_layer().
%   ML [ni, nj]: use a pre-computed mixed layer pressure [dbar] or depth [m]
%   
%   nx [1, 1]: size of ni.
%   ny [1, 1]: size of nj.
%   DXCvec : distance between grid points (i,j) and (i-1,j), always called dx
%   DYCsc : distance between  grid points (i,j) and (i,j-1), always called dy
%   DXGvec : distance between grid points (i+1/2,j) and (i-1/2,j)
%   DYGsc : distance between grid points (i,j+1/2) and (i,j-1/2)
%   RACvec : area of the grid in tracer grid.
%   XCvec : longitude on the T grid.
%   YCvec : latitude on the T grid.
%   DAMP: damping rate of the iterative update of the height of the surface 
%   x0: latitude of the reference cast
%   y0: longitude of the reference cast
%
%
% --- References:

% Lang, Stanley, McDougall: Spurious dianeutral advection and methods for 
%  its minimization.
%
% Stanley, McDougall, Barker 2021: Algorithmic improvements to finding
%  approximately neutral surfaces, Journal of Advances in Earth System
%  Modelling, 13(5).


% Author(s) : Yandong Lang, Geoff Stanley
% Email     : yandong.lang@unsw.edu.au, g.stanley@unsw.edu.au
% Email     : aaronlangyandong@gmail.com, geoffstanley@gmail.com


% --- Notes on the code:
% Upper case letters, e.g. S, denote 3D scalar fields [nk,ni,nj]
% Lower case letters, e.g. s, denote 2D scalar fields    [ni,nj]
% Developmental things are marked with a comment "DEV"


%% Defining functions
im1 = @(F) circshift(F, [+1 0]);
jm1 = @(F) circshift(F, [0 +1]);
ip1 = @(F) circshift(F, [-1 0]);
jp1 = @(F) circshift(F, [0 -1]);
lead1 = @(x) reshape(x, [1 size(x)]);

rms = @(f) sqrt(mean(f(:).^2));
nanrms = @(f) sqrt(nanmean(f(:).^2));

ni = OPTS.nx; nj = OPTS.ny;

%% Grid information
DXC = OPTS.DXCvec + zeros(OPTS.nx, OPTS.ny); %dx
DXG = OPTS.DXGvec + zeros(OPTS.nx, OPTS.ny); %dxT
DYC = OPTS.DYCsc  + zeros(OPTS.nx, OPTS.ny); %dy
DYG = OPTS.DYGsc  + zeros(OPTS.nx, OPTS.ny); %dyT
RAC = OPTS.RACvec + zeros(OPTS.nx, OPTS.ny);
RAX = DXC .* DYG; % area on (i-1/2, j), u cell
RAY = DXG .* DYC; % area on (i, j-1/2), v cell
sqrtAREA = sqrt(RAC);
sqrtAREAX = sqrt(RAX);
sqrtAREAY = sqrt(RAY);

A5 = grid_adjacency([ni,nj], 5, [1; 1]); % this code currently needs doubly periodic domain... easy to fix
A4 = A5(1:4,:);

BotK = max(1,squeeze(sum(isfinite(S))));
Z_bot = Z(BotK);

%% Interpolation
interpfn = OPTS.INTERPFN;
% interpfn = @ppc_linterp;
% interpfn = @ppc_pchip;
SppZ = interpfn(Z, S); % S as a function of Z
TppZ = interpfn(Z, T);
UppZ = interpfn(Z, U); % knowing Z is 1D vector, horizontally uniform
VppZ = interpfn(Z, V); % knowing Z is 1D vector, horizontally uniform

%% Ref cast
if OPTS.data_cube
    i0 = OPTS.i0;
    j0 = OPTS.j0;
    I0 = sub2ind([OPTS.nx,OPTS.ny], i0, j0);
else
    x0 = OPTS.x0;
    y0 = OPTS.y0;
    [~,i0] = min(abs(OPTS.XCvec - x0));
    [~,j0] = min(abs(OPTS.YCvec - y0));
    I0 = sub2ind([OPTS.nx,OPTS.ny], i0, j0);
end

%% Diagnostics
n_iter = OPTS.ITER;

d = struct();
d.ehel_rms = nan(n_iter,1);
d.ehel_s2_rms = nan(n_iter,1);
d.ehelx_rms = nan(n_iter,1);
d.ehely_rms = nan(n_iter,1);
d.dz_rms = nan(n_iter,1);
d.sx_rms = nan(n_iter,1);
d.sy_rms = nan(n_iter,1);
d.s2x_rms = nan(n_iter,1);
d.s2y_rms = nan(n_iter,1);
d.s2_rms = nan(n_iter,1);
d.dz_zd = nan(n_iter, ni, nj);
d.clocktime = nan(n_iter,1);
d.mean_z = nan(n_iter, 1);
d.mean_eos = nan(n_iter, 1);

%% Run omega_hel
SHEAR = OPTS.SHEAR;
damp = OPTS.DAMP;

for iter = 1 : n_iter
  iter_tic = tic;
  
  if iter == 1 || SHEAR % quick way to make 3D velocity have zero shear.
    u = ppc_val_mex(Z, UppZ, (z + im1(z)) / 2, 0);
    v = ppc_val_mex(Z, VppZ, (z + jm1(z)) / 2, 0);
  end
  uz = ppc_val_mex(Z, UppZ, (z + im1(z)) / 2, 1);
  vz = ppc_val_mex(Z, VppZ, (z + jm1(z)) / 2, 1);
  
  Tz = ppc_val_mex(Z, TppZ, z, 1);
  Sz = ppc_val_mex(Z, SppZ, z, 1);
  
  Tzz = ppc_val_mex(Z, TppZ, z, 2);
  Szz = ppc_val_mex(Z, SppZ, z, 2);
  
  % --- Compute slope errors
  tolz = 1e-8;
  [dzi, dzj] = ntp_slope_mex(SppZ, TppZ, Z, z, tolz, 1, 1);
  dzi = -dzi; dzj = -dzj; % undo minus sign in ntp_slope
  sx = -(dzi - (z - im1(z)) ) ./ OPTS.DXCvec; % minus sign to account for z > 0
  sy = -(dzj - (z - jm1(z)) ) ./ OPTS.DYCsc;
  
  mean_z = nanmean(z(:));
  [s,t] = ppc_val2(Z,SppZ,TppZ,lead1(z));
  mean_eos = nanmean(eos(s(:), t(:), z(:)));
  
  %% Matrix
  mode = OPTS.MODE;
  if     mode == 1 % omega_hel
      [dz, esp, N] = omega_hel_matsolve(z, Z, SppZ, TppZ, u, v, uz, vz, sx, sy, dzi, dzj, DXG, DYG, DXC, DYC, RAX, RAY, RAC, sqrtAREA, i0, j0, I0, A4, A5, OPTS);
  elseif mode == 2 % omega_s
      [dz, N]       = omega_s_matsolve(z, Z, SppZ, TppZ, sqrtAREAX, sqrtAREAY, sx, sy, dzi, dzj, i0, j0, I0, A4, OPTS);
  elseif mode == 3 % omega_1.5
      [dz,ehelx,ehely,N] = omega_1p5_matsolve(z, Z, SppZ, TppZ, u, v, uz, vz, sx, sy, dzi, dzj, i0, j0, I0, A4, OPTS);
  elseif mode == 4 % omega_s^2
      [dz, s2, N]        = omega_s2_matsolve2(z, Z, SppZ, TppZ, sx, sy, dzi, dzj, DXG, DYG, RAC, sqrtAREA, i0, j0, I0, A4, A5, OPTS);
  elseif mode == 5 % combining omega_hel and omega_s^2
      [dz, esp, N] = omega_hels2_matsolve(z, Z, SppZ, TppZ, u, v, uz, vz, sx, sy, dzi, dzj, DXG, DYG, RAX, RAY, RAC, sqrtAREA, i0, j0, I0, A4, A5, OPTS);
  elseif mode == 6 % omega_s2xy minimizing sx^2 and sy^2 seperately
      [dz, s2x, s2y, N] = omega_s2xy_matsolve(z, Z, SppZ, TppZ, sqrtAREAX, sqrtAREAY, sx, sy, dzi, dzj, i0, j0, I0, A4, OPTS);
  elseif mode == 7 % omega_hel combines with omega_1.5
      [dz,ehelx,ehely, esp, N] = omega_hel1p5_matsolve(z, Z, SppZ, TppZ, u, v, uz, vz, sx, sy, dzi, dzj, DXG, DYG, RAX, RAY, RAC, sqrtAREA, i0, j0, I0, A4, A5, OPTS);
  elseif mode == 8 % omega_hel combines with omega_s2xy
      [dz,s2x,s2y, esp, N] = omega_hels2xy_matsolve(z, Z, SppZ, TppZ, sqrtAREAX, sqrtAREAY,  u, v, uz, vz, sx, sy, dzi, dzj, DXG, DYG, DXC, DYC, RAX, RAY, RAC, sqrtAREA, i0, j0, I0, A4, A5, OPTS);
  elseif mode == 9 % omega_hel weighted by Tz
      [dz, ehelTz, N] = omega_hel_Tz_matsolve(z, Z, SppZ, TppZ, Tz, Tzz, u, v, uz, vz, sx, sy, dzi, dzj, DXG, DYG, DXC, DYC, RAX, RAY, RAC, sqrtAREA, i0, j0, I0, A4, A5, OPTS);
  elseif mode == 10 % omega_hel weighted by Sz
      [dz, ehelSz, N] = omega_hel_Sz_matsolve(z, Z, SppZ, TppZ, Sz, Szz, u, v, uz, vz, sx, sy, dzi, dzj, DXG, DYG, DXC, DYC, RAX, RAY, RAC, sqrtAREA, i0, j0, I0, A4, A5, OPTS);
  end
    
  % --- Update the surface
  z = z + damp * dz;
  
  z(z < Z(1) | z > Z_bot) = nan;

  %% --- Diagnostics
  if mode == 1
      d.ehel_rms(iter) = rms(esp(:));
  elseif mode == 5
      d.ehel_rms(iter) = rms(esp(:));
  elseif mode == 3
      d.ehelx_rms(iter) = nanrms(ehelx(:));
      d.ehely_rms(iter) = nanrms(ehely(:)); 
  elseif mode == 4
      d.s2_rms(iter) = nanrms(s2(:));
  elseif mode == 6
      d.s2x_rms(iter) = nanrms(s2x(:));
      d.s2y_rms(iter) = nanrms(s2y(:));
  elseif mode == 7
      d.ehel_rms(iter) = root_mean_square(esp(:));
      d.ehelx_rms(iter) = nanrms(ehelx(:));
      d.ehely_rms(iter) = nanrms(ehely(:));
  elseif mode == 8
      d.ehel_s2_rms(iter) = root_mean_square(esp(:)) + nanrms(s2x(:))+ nanrms(s2y(:));
  elseif mode == 9
      d.ehelTz_rms(iter) = root_mean_square(ehelTz(:));
  elseif mode == 10
      d.ehelSz_rms(iter) = root_mean_square(ehelSz(:));
  end
  d.dz_rms(iter) = nanrms(abs(dz(:)));
  d.sx_rms(iter) = nanrms(sx(:));
  d.sy_rms(iter) = nanrms(sy(:));
  d.num_casts(iter) = N;
  d.dz_2d(iter,:,:) = dz;
  d.mean_z(iter) = mean_z;
  d.mean_eos(iter) = mean_eos;
  d.clocktime(iter) = toc(iter_tic);
  
  % diagnostic sentence
  if mode == 2
      fprintf('Iter %2d: (%6.2f sec), log10(|sx|_2) = %.6f, log10(|sy|_2) = %.6f, log10(|Δz|_2) = %.2f, mean depth = %.2f, mean(eos) = %.6e, # casts = %4d\n', ...
          iter, d.clocktime(iter), log10(d.sx_rms(iter)), log10(d.sy_rms(iter)), log10(d.dz_rms(iter)), d.mean_z(iter), d.mean_eos(iter), N);
  elseif mode == 1
      fprintf('Iter %2d: (%6.2f sec), log10(|ehel|_2) = %.6f, log10(|Δz|_2) = %.2f, mean depth = %.2f, mean(eos) = %.6e, # casts = %4d\n', ...
          iter, d.clocktime(iter), log10(d.ehel_rms(iter)), log10(d.dz_rms(iter)), d.mean_z(iter), d.mean_eos(iter), N);
  elseif mode == 5
      fprintf('Iter %2d: (%6.2f sec), log10(|ehel|_2) = %.6f, log10(|Δz|_2) = %.2f, mean depth = %.2f, mean(eos) = %.6e, # casts = %4d\n', ...
          iter, d.clocktime(iter), log10(d.ehel_rms(iter)), log10(d.dz_rms(iter)), d.mean_z(iter), d.mean_eos(iter), N);
  elseif mode == 3
      fprintf('Iter %2d: (%6.2f sec), log10(|ehelx|_2) = %.6f, log10(|ehely|_2) = %.6f, log10(|Δz|_2) = %.2f, mean depth = %.2f, mean(eos) = %.6e, # casts = %4d\n', ...
          iter, d.clocktime(iter), log10(d.ehelx_rms(iter)), log10(d.ehely_rms(iter)), d.mean_z(iter), d.mean_eos(iter), log10(d.dz_rms(iter)), N);
  elseif mode == 4
      fprintf('Iter %2d: (%6.2f sec), log10(|s^2|_2) = %.6f, log10(|Δz|_2) = %.2f, mean depth = %.2f, mean(eos) = %.6e, # casts = %4d\n', ...
          iter, d.clocktime(iter), log10(d.s2_rms(iter)), log10(d.dz_rms(iter)), d.mean_z(iter), d.mean_eos(iter), N); 
  elseif mode == 6
      fprintf('Iter %2d: (%6.2f sec), log10(|s2_x|_2) = %.6f, log10(|s2y|_2) = %.6f, log10(|Δz|_2) = %.2f, mean depth = %.2f, mean(eos) = %.6e, # casts = %4d\n', ...
          iter, d.clocktime(iter), log10(d.s2x_rms(iter)), log10(d.s2y_rms(iter)), log10(d.dz_rms(iter)), d.mean_z(iter), d.mean_eos(iter), N);
  elseif mode == 7
      fprintf('Iter %2d: (%6.2f sec), log10(|ehel|_2) = %.6f, log10(|ehelx|_2) = %.6f, log10(|ehely|_2) = %.6f, log10(|Δz|_2) = %.2f, mean depth = %.2f, mean(eos) = %.6e, # casts = %4d\n', ...
          iter, d.clocktime(iter), log10(d.ehel_rms(iter)), log10(d.ehelx_rms(iter)), log10(d.ehely_rms(iter)), log10(d.dz_rms(iter)), d.mean_z(iter), d.mean_eos(iter), N);
  elseif mode == 8
      fprintf('Iter %2d: (%6.2f sec), log10(|ehel+s2|_2) = %.6f, log10(|Δz|_2) = %.2f, mean depth = %.2f, mean(eos) = %.6e, # casts = %4d\n', ...
          iter, d.clocktime(iter), log10(d.ehel_s2_rms(iter)), log10(d.dz_rms(iter)), d.mean_z(iter), d.mean_eos(iter), N);
  elseif mode == 9
      fprintf('Iter %2d: (%6.2f sec), log10(|ehelTz|_2) = %.6f, log10(|Δz|_2) = %.2f, mean depth = %.2f, mean(eos) = %.6e, # casts = %4d\n', ...
          iter, d.clocktime(iter), log10(d.ehelTz_rms(iter)), log10(d.dz_rms(iter)), d.mean_z(iter), d.mean_eos(iter), N);
  elseif mode == 10
      fprintf('Iter %2d: (%6.2f sec), log10(|ehelSz|_2) = %.6f, log10(|Δz|_2) = %.2f, mean depth = %.2f, mean(eos) = %.6e, # casts = %4d\n', ...
          iter, d.clocktime(iter), log10(d.ehelSz_rms(iter)), log10(d.dz_rms(iter)), d.mean_z(iter), d.mean_eos(iter), N);
  end

%% --- Show Figures
FIGS_SHOW = OPTS.FIGS_SHOW;
  if FIGS_SHOW
    if mod(iter,12) == 1
      hf = figure('Position',[20, 20, 1000, 800], 'Name', 'phi''');
    end
    ax = subplot(4, 3, mod(iter-1,12)+1, 'Parent', hf );
    if nj > ni % p, and phi, are probably lon p lat
      %imagesc(ax, phi)
      imagesc(ax, dz)
    else       % p, and phi, are probably lat p lon
      %imagesc(ax, phi.')
      imagesc(ax, dz')
    end
    ax.YDir = 'normal';
    shading('flat')
    title(sprintf('%d: log10(|Δz|_2) = %.2f', iter, log10(d.dz_rms(iter))), 'fontsize',10);
    caxis(prctile(dz(:), [1 99]));
    colorbar();
    drawnow()
  end
  
end

zns = z;
[sns,tns] = ppc_val2(Z,SppZ,TppZ,lead1(zns));

% Diagnostics for final surface
ehelx = u .* sx .* (OPTS.DXCvec .* OPTS.DYGsc);
ehely = v .* sy .* (OPTS.DXGvec .* OPTS.DYCsc);
ehelx(isnan(ehelx)) = 0;
ehely(isnan(ehely)) = 0;
esp = (ehelx + ip1(ehelx) + ehely + jp1(ehely)) ./ (2 * OPTS.RACvec);
% d.z_map(:,:,iter+1) = z;
% d.ehel_map(:,:,iter+1) = ehel;
end