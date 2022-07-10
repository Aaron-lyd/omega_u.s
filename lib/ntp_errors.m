function [ex,ey,sx,sy] = ntp_errors(s, t, p, dx, dy, use_s_t, centre, wrap, eoses, grav, S, T, P)
%NTP_ERRORS  The epsilon neutrality error and slope error from the neutral tangent plane
%
%
% [ex,ey] = ntp_errors(s,t,p,dx,dy,use_s_t,centre,wrap,eoses)
% computes the zonal and meridional neutrality errors, ex and ey, for an
% approximately neutral surface on which the practical / Absolute salinity
% is s, the potential / Conservative temperature is t, and the pressure or
% depth is p.
% The zonal and meridional distances between neighbouring grid points are
% dx and dy.
% If use_s_t is true, errors are calculated using s and t differences
% multiplied by the partial derivatives of in-situ density w.r.t s and t,
% respectively; if use_s_t is false, errors are calculated using in-situ
% density differences and p differences multiplied by the partial
% derivative of in-situ density w.r.t. p.
% If centre is true, errors are calculated using centred differences, and
% results are on the grid points of s, t and p; otherwise, backward
% differences are used and results for ex and ey are located zonally and
% meridionally, respectively, midway between s, t, and p grid points.
% Data is treated periodic in the i'th (i=1 or 2) dimension of p if and
% only if wrap(i) is true.
%
% If eoses is empty, then the equation of state for the in-situ density is
% given by eos.m in the path.  If use_s_t is true, eos_s_t.m must also
% exist in the path and return the partial derivatives of in-situ density
% w.r.t s and t; otherwise, eos_p.m must also exist in the path and return
% the partial derivative of in-situ density w.r.t p. The inputs to these
% eos*.m functions must be (s, t, p).  Note, eos*.m can involve specific
% volume instead of in-situ density, which merely changes the units of
% [ex,ey].
%
% Otherwise, eoses must be a 2 element cell array that allows overriding
% eos.m and eos_s_t.m / eos_p.m on MATLAB's path. eoses{1} is a function
% handle to the equation of state, and eoses{2} is a function handle to the
% necessary partial derivative(s) of the equation of state in its second
% entry.
%
% [ex,ey,sx,sy] = ntp_errors(s,t,p,dx,dy,use_s_t,centre,wrap,eoses,grav,S,T,P)
% also computes the zonal and meridional slope errors, sx and sy, for an
% ocean of practical / Absolute salinity S and potential / Conservative
% temperature T at datasites where the pressure or depth is P.
% Provide {} for eoses if the default files on MATLAB's path are desired.
%
%
% For a non-Boussinesq ocean, P and p are pressure [dbar].  When sx, sy are
% requested, the gravitational acceleration grav must be given.
%
% For a Boussinesq ocean, P and p are depth [m, positive].  When sx, sy are
% requested, grav must be the empty vector [].
%
%
% --- Input:
% s     [ni, nj]: practical / Absolute salinity on the surface
% t     [ni, nj]: potential / Conservative temperature on the surface
% p     [ni, nj]: pressure [dbar] or depth [m, positive] of the surface
% dx    [ni, nj]: Zonal      grid distance [m] between S(:,i,j) and S(:,i-1,j)
% dy    [ni, nj]: Meridional grid distance [m] between S(:,i,j) and S(:,i,j-1)
% use_s_t [1, 1]: true to compute ex,ey using s and t differences,
%                 false to use eos(s,t,p) and p differences.
% centre  [1, 1]: true to compute outputs on the central grid,
%                 false on the half grids
% wrap    [2, 1]: wrap(i) is true when the data is periodic in i'th
%                 dimension of p
% eoses   {2, 1}: cell array of function handles to eos.m and eos_s_t.m (if
%                 use_s_t == true) or eos_p.m (if use_s_t == false).
% grav    [1, 1]: gravitational acceleration [m s^-2]
% S [nk, ni, nj]: practical / Absolute salinity
% T [nk, ni, nj]: potential / Conservative temperature
% P [nk, ni, nj] or [nk, 1]: pressure [dbar] or depth [m, positive]
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in longitude,
%       nj is the number of data points in latitude.
%
% Note: data can be arranged as latitude by longitude.  Simply swap dx and
%       dy, adjust wrap as needed, and note the outputs ex and ey will also
%       be swapped, as will be sx and sy.
%
% Note: physical units for s,t,p,S,T,P are determined by eos*.m.
%
% Note: P must be increasing along its first dimension.
%
% Note: dx and dy can have any dimension replaced by a singleton dimension.
%
%
% --- Output:
% ex [ni, nj]: zonal epsilon neutrality error [(kg m^-3) / m or (m^3 / kg) / m]
% ey [ni, nj]: meridional epsilon neutrality error [(kg m^-3) / m or (m^3 / kg) / m]
% sx [ni, nj]: zonal slope error [dimensionless]
% sy [ni, nj]: meridional slope error [dimensionless]
%
%
% --- Mathematical form of ex, ey, sx, sy:
% The neutrality error (usually denoted by epsilon) is here taken as one of
% two mathematically equivalent forms:
%   [ex,ey] = rs grad s + rt grad t
% when use_s_t is true, or
%   [ex,ey] = grad r - rx grad p
% when use_s_t is false.
% Here, the pressure or depth of the surface is p; s and t are the values
% of S and T on the surface projected onto a sphere; r = eos(s,t,p) is the
% in-situ density on the surface; rs, rt, and rx are the partial
% derivatives of in-situ density with respect to S, T, and P respectively,
% evaluated on the surface. Note that grad s is the same as grad_a S
% evaluated on the surface, where grad_a is the projected non-orthogonal
% gradient. See discussion of the under-tilde operator in Stanley (2019).
%
% The above definition is equivalent to Stanley (2019) Eq. (28), but
% differs from the standard definition by a factor of r, the in-situ
% density. The standard definition (McDougall and Jackett 1988, p. 162) is
% [ex,ey] / r.
%
% The slope error is related to the neutrality error by
%   [sx,sy] = g / (r N2) [ex,ey] = -1 / (d sigma / d z) [ex,ey]
% where g is the gravitational acceleration and N2 is the square of the
% vertical stratification, and sigma is the locally referenced potential
% density. This differs from the standard relation (see McDougall and
% Jackett 1988 Eq. (14), or Klocker and McDougall 2009, Eq. (10)) by a
% factor of r because of our definition for [ex,ey] above.
%
%
% --- References:
% Klocker, A., McDougall, T. J. & Jackett, D. R. A new method for forming
% approximately neutral surfaces. Ocean Science 5, 155?172 (2009).
%
% McDougall, T. J. & Jackett, D. R. On the helical nature of neutral
% trajectories in the ocean. Progress in Oceanography 20, 153?183 (1988).
%
% Stanley, G.J., 2019. Neutral surface topology. Ocean Modelling 138,
% 88–106. https://doi.org/10.1016/j.ocemod.2019.01.008

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


% --- Input checks, set parameters
narginchk(8,13)

[nx,ny] = size(s);
is2D = @(F) ismatrix(F) && all(size(F) == [nx,ny]);
is2Ds = @(F) ismatrix(F) && ismember(size(F,1), [1,nx]) && ismember(size(F,2), [1,ny]);
is3Ds = @(F) ismember(size(F,2), [1,nx]) && ismember(size(F,3), [1,ny]);
assert(is2D(t)  , 'the dimensions of t must match those of s');
assert(is2D(p)  , 'the dimensions of p must match those of s');
assert(is2Ds(dx), 'the dimensions of dx must match those of s (or be singletons)');
assert(is2Ds(dy), 'the dimensions of dy must match those of s (or be singletons)');
assert(isscalar(use_s_t), 'use_s_t must be a (logical) scalar');
assert(isscalar(centre), 'centre must be a (logical) scalar');
assert(numel(wrap) == 2, 'wrap must be a 2 element (logical) vector');


if nargout > 2
  assert(nargin == 13, 'If sx, sy are requested, must provide grav, S, T, P');
  nonBOUSSINESQ = ~isempty(grav); % grav, needed in hydrostatic balance when non-Boussinesq, is given.
  if nonBOUSSINESQ
    assert(isscalar(grav), 'grav, if non-empty, must be a scalar');
  end
  assert(is3Ds(S), 'The last two dimensions of S must match those of s (or be singletons).')
  assert(is3Ds(T), 'The last two dimensions of T must match those of s (or be singletons).')
  assert(is3Ds(P), 'The last two dimensions of P must match those of s (or be singletons).')
  assert(size(T,1) == size(S,1) && size(S,1) == size(P,1), 'The first dimensions of S and T and P must match.');
end

if nargin < 9 || isempty(eoses)
  % Get eos functions from MATLAB path
  assert(~isempty(which('eos')), 'Cannot find eos.m in the path or current directory.');
  eos = @eos;
  if use_s_t
    assert(~isempty(which('eos_s_t')), 'Cannot find eos_s_t.m in the path or current directory.');
    eos_s_t = @eos_s_t;
  else
    assert(~isempty(which('eos_p')), 'Cannot find eos_p.m in the path or current directory.');
    eos_p = @eos_p;
  end
else
  % Load eos functions from input cell array
  assert(iscell(eoses) && numel(eoses) == 2, 'eoses, if provided, must be a cell array of length 2');
  assert(isa(eoses{1},'function_handle'), 'eoses{1} must be a function handle.')
  assert(isa(eoses{2},'function_handle'), 'eoses{2} must be a function handle.')
  eos = eoses{1};
  if use_s_t
    eos_s_t = eoses{2};
  else
    eos_p = eoses{2};
  end
end

eos_is_dens = (eos(34.5, 3, 1000) > 1);
lead1 = @(p) reshape(p, [1 size(p)]);
Pa2db = 1e-4;

if centre
  % Get errors by centred differences on the T grid
  dx = dx * 2;
  dy = dy * 2;
  A_X = @(F) (im1(F) + ip1(F)) / 2;
  D_X = @(F)  ip1(F) - im1(F);
  A_Y = @(F) (jm1(F) + jp1(F)) / 2;
  D_Y = @(F)  jp1(F) - jm1(F);
else
  % Get errors by backward differences, and leave on the U, V grid.
  A_X = @(F) (F + im1(F)) / 2;
  D_X = @(F)  F - im1(F);
  A_Y = @(F) (F + jm1(F)) / 2;
  D_Y = @(F)  F - jm1(F);
end

if ~use_s_t || nargout >= 3
  r = eos(s, t, p); % in-situ density (or specific volume)
end

% Begin calculations for errors in x direction
sa = A_X(s);
ta = A_X(t);
pa = A_X(p);
if use_s_t
  [rsa,rta] = eos_s_t(sa,ta,pa);
  ex = (rsa .* D_X(s) + rta .* D_X(t)) ./ dx;
else
  rpa = eos_p(sa, ta, pa);
  ex = (D_X(r) - rpa .* D_X(p)) ./ dx;
end
% Note, the s & t vs r & p forms of ex above are virtually identical: this
% equivalence is accurate to third order. In particular, note that for
% centred differences, rsa, rta, and rpa need to be evaluated at (sa, ta,
% xa) rather than at the central grid point (s, t, p).
if nargout < 2; return; end

% Begin calculations for errors in y direction.
sa = A_Y(s);
ta = A_Y(t);
pa = A_Y(p);
if use_s_t
  [rsa,rta] = eos_s_t(sa,ta,pa);
  ey = (rsa .* D_Y(s) + rta .* D_Y(t)) ./ dy;
else
  rpa = eos_p(sa, ta, pa);
  ey = (D_Y(r) - rpa .* D_Y(p)) ./ dy;
end
if nargout < 3; return; end

% Start working on slope errors

% Set threshold for minimum magnitude of d(sigma)/dz, where sigma =
% potential density or potential specific volume.
% N^2 = -(g / rho) * d(pot dens)/dz = +(g / v) * d(pot spec vol)/dz,   where v = 1/rho
N_min = 2e-4;
thresh = N_min^2 * eos(34.5, 3, 1000) / 9.81; % > 0.  Use canonical Earth gravity for this scaling.


% Calculate locally referenced potential density (or specific volume)
sigma = eos(S, T, lead1(p));

% Take derivative of sigma w.r.t pressure or depth (as a SPATIAL coordinate)
[~, dsigmadz] = pchipdqn(lead1(p), P, sigma);
if nonBOUSSINESQ
  % P is pressure, not depth.
  % Convert d(sigma)/dp into d(sigma)/d|z| using hydrostatic balance
  if eos_is_dens
    dsigmadz = dsigmadz .* (Pa2db * grav * r);
  else
    dsigmadz = dsigmadz .* (Pa2db * grav ./ r);
  end
end
% Now dsigmadz is d(sigma)/d|z| whether P is pressure or inverted depth
% N.B. don't multiply by -1, because P increases downwards even if it is depth.

% Apply threshold:
if eos_is_dens
  dsigmadz(dsigmadz < thresh) = thresh;
else
  dsigmadz(dsigmadz > -thresh) = -thresh;
end

% Slope error in p direction
sx = ex ./ A_X(dsigmadz);


% Slope error in y direction
sy = ey ./ A_Y(dsigmadz);


  function out = im1(in)
    out = circshift(in, [+1 0]);
    if ~wrap(1)
      out(1,:) = nan;
    end
  end

  function out = jm1(in)
    out = circshift(in, [0 +1]);
    if ~wrap(2)
      out(:,1) = nan;
    end
  end

  function out = ip1(in)
    out = circshift(in, [-1 0]);
    if ~wrap(1)
      out(end,:) = nan;
    end
  end

  function out = jp1(in)
    out = circshift(in, [0 -1]);
    if ~wrap(2)
      out(:,end) = nan;
    end
  end

end
