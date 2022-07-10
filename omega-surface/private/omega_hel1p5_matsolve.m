function [dz,ehelx,ehely, ehel, N] = omega_hel1p5_matsolve(z, Z, SppZ, TppZ, u, v, uz, vz, sx, sy, dzi, dzj, DXG, DYG, RAX, RAY, RAC, sqrtAREA, i0, j0, I0, A4, A5, OPTS)

%%
STRAT = OPTS.STRAT;
SHEAR = OPTS.SHEAR;
STRAT_MIN = 0;
LM = OPTS.LM;
flat = @(F) F(:);
im1 = @(F) circshift(F, [+1 0]);
jm1 = @(F) circshift(F, [0 +1]);
ip1 = @(F) circshift(F, [-1 0]);
jp1 = @(F) circshift(F, [0 -1]);
[ni,nj] = size(z);

u(isnan(u)) = 0; % DEV.  Will be handled by goodx having a u ~= 0 condition, later.
v(isnan(v)) = 0; % DEV.  Will be handled by goody having a v ~= 0 condition, later.

%% Stratification
if STRAT
    % Calculate stratification terms
    z_mn = (z + im1(z)) / 2;
    [sm, tm] = ppc_val2_mex(Z, SppZ, TppZ, z_mn + dzi/2, 0);
    [sn, tn] = ppc_val2_mex(Z, SppZ, TppZ, ip1(z_mn - dzi/2), 0);
    sn = im1(sn); tn = im1(tn);
    [szm, tzm] = ppc_val2_mex(Z, SppZ, TppZ, z_mn + dzi/2, 1);
    [szn, tzn] = ppc_val2_mex(Z, SppZ, TppZ, ip1(z_mn - dzi/2), 1);
    szn = im1(szn); tzn = im1(tzn);
    
    [rsm, rtm] = eos_s_t(sm, tm, z_mn);
    rzm        = eos_x(  sm, tm, z_mn);
    [rsn, rtn] = eos_s_t(sn, tn, z_mn);
    rzn        = eos_x(  sn, tn, z_mn);
    
    H_imh = ((rsm .* szm - rsn .* szn) ...
        +      (rtm .* tzm - rtn .* tzn) ...
        +      (rzm        - rzn       )) ...
        ./    ((rsm .* szm + rsn .* szn) ...
        +      (rtm .* tzm + rtn .* tzn) + STRAT_MIN);
    
    % Repeat for 2nd dimension
    z_mn = (z + jm1(z)) / 2;
    [sm, tm] = ppc_val2_mex(Z, SppZ, TppZ, z_mn + dzj/2, 0);
    [sn, tn] = ppc_val2_mex(Z, SppZ, TppZ, jp1(z_mn - dzj/2), 0);
    sn = jm1(sn); tn = jm1(tn);
    [szm, tzm] = ppc_val2_mex(Z, SppZ, TppZ, z_mn + dzj/2, 1);
    [szn, tzn] = ppc_val2_mex(Z, SppZ, TppZ, jp1(z_mn - dzj/2), 1);
    szn = jm1(szn); tzn = jm1(tzn);
    
    [rsm, rtm] = eos_s_t(sm, tm, z_mn);
    rzm        = eos_x(  sm, tm, z_mn);
    [rsn, rtn] = eos_s_t(sn, tn, z_mn);
    rzn        = eos_x(  sn, tn, z_mn);
    
    H_jmh = ((rsm .* szm - rsn .* szn) ...
        +      (rtm .* tzm - rtn .* tzn) ...
        +      (rzm        - rzn       )) ...
        ./    ((rsm .* szm + rsn .* szn) ...
        +      (rtm .* tzm + rtn .* tzn) + STRAT_MIN);
end

%% --- Compute u . s.  See emails from Geoff to Aaron on July 16 and 19, 2021.
%u . s = 1 / (2Ai,j) ( ui-1/2,j sxi-1/2,j Δxi-1/2,j Δyi-1/2,j + ui+1/2,j sxi+1/2,j Δxi+1/2,j Δyi+1/2,j + vi,j-1/2 syi,j-1/2 Δyi,j-1/2 Δxi,j-1/2 + vi,j+1/2 syi,j+1/2 Δyi,j+1/2 Δxi,j+1/2 )
ehelx = u .* sx .* (OPTS.DXCvec .* OPTS.DYGsc);
ehely = v .* sy .* (OPTS.DXGvec .* OPTS.DYCsc);
ehelx(isnan(ehelx)) = 0; % also ehelx == 0  when u == 0, so ehelx == 0 when ~goodx
ehely(isnan(ehely)) = 0; % also ehely == 0  when v == 0, so ehely == 0 when ~goody
ehel = (ehelx + ip1(ehelx) + ehely + jp1(ehely)) ./ (2 * OPTS.RACvec);

%% Find casts where at least one neighbour provides a valid connection
goodx = isfinite(sx) & u ~= 0;
goody = isfinite(sy) & v ~= 0;
goodx_ip1 = ip1(goodx);
goody_jp1 = jp1(goody);
% valid = isfinite(z) & (goodx ~= 0 | goody ~= 0 | goodx_ip1 ~= 0 | goody_jp1 ~= 0);
% [qu, qt] = bfs_conncomp1(valid, A4, I0);

% New method to find casts where at least one neighbour provides a valid connection
A4_ = A4;
wallval = ni*nj+1;
A4_(1,~goody) = wallval;
A4_(2,~goodx) = wallval;
A4_(3,~ip1(goodx)) = wallval;
A4_(4,~jp1(goody)) = wallval;
[qu, qt] = bfs_conncomp1(isfinite(z), A4_, I0);

%% Matrix omega_hel
m = sort(qu(1:qt));
N = length(m);  % number of water columns with at least one valid equation

remap = zeros(ni, nj);  % map from 2D space to 1D vector of water columns
remap(m) = 1 : N;

% -+-----> j
%  | . 2 .
%  | 1 5 4
%  v . 3 .
%  i
r = reshape(repelem(1 : N, 5), 5, N);
c = remap(A5(:,m));

% Prepare array of valid connections from the central cast to adjacent
% casts. A connection is valid when the NTP between the pair was
% successfully found, i.e. is not NaN.  (Note here that the fifth row of
% this definition is inconsequential.)  This can be slightly different
% from the definition of a valid pair being when the surface is wet (not
% in/outcropped) at both casts, since the NTP can be steeper than the
% surface.  If we defined valid connections by the second criteria (wet
% surface casts), we could just code this as: >> good = c > 0;

% Select only pairs of casts where the NTP exists and the velocity is non-zero.
% For a realistic velocity field, the velocity is zero only when the surface
% average depth between two casts is deeper than the topography between
% those two casts, so the velocity is interpolated to a depth at which there
% is a vertical wall, where the normal velocity is zero.  In this case, a
% small adjustment up or down will still yield zero velocity, so zero
% contribution to ehel, so we should ignore those contributions.
good = [goody(m)'; goodx(m)'; goodx_ip1(m)'; goody_jp1(m)'; true(1,N)];
assert(all(any(good(1:4,:))), 'Found cast with all four equations invalid')

% Do (u,v) . grad(phi) terms
im = u .* DYG; ip = ip1(im);
jm = v .* DXG; jp = jp1(jm);
% val = [-jm(m)'; -im(m)'; ip(m)'; jp(m)'; zeros(1, N)] ./ (2 * RAC(m)');
val = [-jm(m)'; -im(m)'; ip(m)'; jp(m)'; zeros(1, N)] ./ (2 * sqrtAREA(m)');
val(~good) = 0; % Ignore invalid connections
val(5,:) = -sum(val(1:4,:), 1); % Note - sign

assert(sum(val(good) == 0) == 0, 'found %d zeros in val(good)', sum(val(good) == 0))
if SHEAR
    % Do (u',v') . s terms
    im = uz .* sx .* RAX; ip = ip1(im);
    jm = vz .* sy .* RAY; jp = jp1(jm);
%     val2 = [jm(m)'; im(m)'; ip(m)'; jp(m)'; zeros(1,N)] ./ (2 * RAC(m)');
    val2 = [jm(m)'; im(m)'; ip(m)'; jp(m)'; zeros(1,N)] ./ (4 * sqrtAREA(m)');
    val2(~good) = 0; % Ignore invalid connections
    val2(5,:) = sum(val2(1:4,:), 1); % note + sign
    val = val + val2;
end

if STRAT
    % Do (u,v) . H terms
    im = u .* H_imh .* DYG; ip = ip1(im);
    jm = v .* H_jmh .* DXG; jp = jp1(jm);
%     val2 = [jm(m)'; im(m)'; ip(m)'; jp(m)'; zeros(1,N)] ./ (2 * RAC(m)');
    val2 = [jm(m)'; im(m)'; ip(m)'; jp(m)'; zeros(1,N)] ./ (2 * sqrtAREA(m)');
    val2(~good) = 0; % Ignore invalid connections
    val2(5,:) = sum(val2(1:4,:), 1); % note + sign
    val = val + val2;
end

mat_hel = sparse(r(good), c(good), val(good), N, N);

rhs_hel = zeros(N,1);
rhs_hel(1:N) = -ehel(m) .* sqrtAREA(m);

% final row values for pinning
rf = rms(val(good));

%% Matrix omega_1.5
good_1p5 = [goodx(m); goody(m); true(1,N)'; true(1,N)'];

DXC = OPTS.DXCvec + zeros(ni, nj);
DYC = OPTS.DYCsc  + zeros(ni, nj);

DXCm = DXC(m);
DYCm = DYC(m);

% Build the RHS of the matrix problem
sxm = sx(m);
sym = sy(m);
um = u(m);
vm = v(m);

ehelx = sxm .* um;
ehely = sym .* vm;

rhs_1p5 = -[ehelx; ehely;0];
neqx = length(ehelx);
neqy = length(ehely);

neq = neqx + neqy; % Number of equations, excluding density conserving equation.

% Build indices for the rows of the sparse matrix
row = [(1:neq)';(1:neq)'];

% Build indices for the columns of the sparse matrix of slope error
c = remap(A4(:,m).'); 
c_0 = flat(1:N); % z'
c_1 = c(:,1); % jm1
c_2 = c(:,2); % im1

col = [c_2; c_1; c_0; c_0];

% Build the values of the sparse matrix
p1x = um;
p1y = vm;

val = [-p1x ./ DXCm; -p1y ./ DYCm; p1x ./ DXCm; p1y ./ DYCm];
    
val(~good_1p5) = 0; % Ignore invalid connections

goodxm = goodx(m);
goodym = goody(m);
if STRAT
    % H terms
    H_imh_m = H_imh(m);
    H_jmh_m = H_jmh(m);
    
    im = p1x .* H_imh_m ./ DXCm;
    im(~goodxm) = 0; % Ignore invalid connections
    ip = im;
    jm = p1y .* H_jmh_m ./ DYCm;
    jm(~goodym) = 0; % Ignore invalid connections
    jp = jm;
    
    val2 = [im; jm; ip; jp];
    val = val + val2;
end

if SHEAR
    uzm = uz(m);
    vzm = vz(m);
    im = uzm .* sxm/2;
    im(~goodxm) = 0; % Ignore invalid connections
    ip = im;
    jm = vzm .* sym/2;
    jm(~goodym) = 0; % Ignore invalid connections
    jp = jm;
    
    val2 = [im; jm; ip; jp];
    val = val + val2;
end

% Build the sparse matrix, with neq rows and N columns
mat_1p5 = sparse(row(good_1p5), col(good_1p5), val(good_1p5), neq, N );

%% Combining matrix

% reference cast
mr = remap(i0,j0);
% weight of ehelu and ehelv
ww = OPTS.UVW;
mat = [mat_hel; ww*mat_1p5];
rhs = [rhs_hel; ww*rhs_1p5];

mat(N+neq+1,mr) = rf;

% mat = [mat_1p5];
% rhs = [rhs_1p5];
% mat(neq+1,mr) = rf;

% Normal Equations
rhs = mat' * rhs;
mat = mat' * mat;

%% Levenberg-Marquardt
if LM > 0
    mat_rms = rms(mat(mat ~= 0));
    assert(isfinite(mat_rms), 'matrix has some NaN''s; cannot proceed.');
    mat = mat + (LM * mat_rms) * speye(N);
end

%% TK
TK = OPTS.TK;
if TK > 0
    % using a high-pass filter, specifically a difference operator: grad phi
    num_eq = sum(goodx(m)) + sum(goody(m));
    r = reshape(repelem(1 : num_eq, 2), 2, num_eq);
    c = [remap(A4(1,m)), remap(A4(2,m)); 1:N, 1:N];
    c = c(:,[goody(m), goodx(m)]);
    
    % uniform grid of unity
    val = repelem([-1; 1], 1, num_eq);
    tik = sparse(r, c, val, num_eq, N);
    
    % Apply regularization.  Note tik' * tik is a Laplacian.  Could build that directly...
    mat = mat + TK * (tik' * tik);
    % Note that adding Tikhonov makes pinning not quite perfect.
    % In one test, the solution drifted 6.36cm at ref cast...  current solution is PIN_RIGID
end

%% solution
sol = mat \ rhs;

dz = zeros(ni, nj);
dz(m) = sol;

end