function [dz, s2x, s2y, N] = omega_s2xy_matsolve(z, Z, SppZ, TppZ, sqrtAREAX, sqrtAREAY, sx, sy, dzi, dzj, i0, j0, I0, A4, OPTS)

%
STRAT = OPTS.STRAT;
H_SIM = OPTS.H_SIM;
STRAT_MIN = 0;
LM = OPTS.LM;
im1 = @(F) circshift(F, [+1 0]);
jm1 = @(F) circshift(F, [0 +1]);
ip1 = @(F) circshift(F, [-1 0]);
jp1 = @(F) circshift(F, [0 -1]);
[ni,nj] = size(z);
flat = @(F) F(:);

LM = OPTS.LM;
TK = OPTS.TK;

%% Stratification terms
% if STRAT
%     % Calculate stratification terms
%     z_mn = (z + im1(z)) / 2;
%     [sm, tm] = ppc_val2_mex(Z, SppZ, TppZ, z_mn + dzi/2, 0);
%     [sn, tn] = ppc_val2_mex(Z, SppZ, TppZ, ip1(z_mn - dzi/2), 0);
%     sn = im1(sn); tn = im1(tn);
%     [szm, tzm] = ppc_val2_mex(Z, SppZ, TppZ, z_mn + dzi/2, 1);
%     [szn, tzn] = ppc_val2_mex(Z, SppZ, TppZ, ip1(z_mn - dzi/2), 1);
%     szn = im1(szn); tzn = im1(tzn);
%     
%     [rsm, rtm] = eos_s_t(sm, tm, z_mn);
%     rzm        = eos_x(  sm, tm, z_mn);
%     [rsn, rtn] = eos_s_t(sn, tn, z_mn);
%     rzn        = eos_x(  sn, tn, z_mn);
%     
%     H_imh = ((rsm .* szm - rsn .* szn) ...
%         +      (rtm .* tzm - rtn .* tzn) ...
%         +      (rzm        - rzn       )) ...
%         ./    ((rsm .* szm + rsn .* szn) ...
%         +      (rtm .* tzm + rtn .* tzn) + STRAT_MIN);
%     
%     % Repeat for 2nd dimension
%     z_mn = (z + jm1(z)) / 2;
%     [sm, tm] = ppc_val2_mex(Z, SppZ, TppZ, z_mn + dzj/2, 0);
%     [sn, tn] = ppc_val2_mex(Z, SppZ, TppZ, jp1(z_mn - dzj/2), 0);
%     sn = jm1(sn); tn = jm1(tn);
%     [szm, tzm] = ppc_val2_mex(Z, SppZ, TppZ, z_mn + dzj/2, 1);
%     [szn, tzn] = ppc_val2_mex(Z, SppZ, TppZ, jp1(z_mn - dzj/2), 1);
%     szn = jm1(szn); tzn = jm1(tzn);
%     
%     [rsm, rtm] = eos_s_t(sm, tm, z_mn);
%     rzm        = eos_x(  sm, tm, z_mn);
%     [rsn, rtn] = eos_s_t(sn, tn, z_mn);
%     rzn        = eos_x(  sn, tn, z_mn);
%     
%     H_jmh = ((rsm .* szm - rsn .* szn) ...
%         +      (rtm .* tzm - rtn .* tzn) ...
%         +      (rzm        - rzn       )) ...
%         ./    ((rsm .* szm + rsn .* szn) ...
%         +      (rtm .* tzm + rtn .* tzn) + STRAT_MIN);
% end

if STRAT
    if H_SIM
        % Calculate stratification terms
        [sm, tm] = ppc_val2_mex(Z, SppZ, TppZ, z, 0);
        sn = im1(sm); tn = im1(tm);
        [szm, tzm] = ppc_val2_mex(Z, SppZ, TppZ, z, 1);
        szn = im1(szm); tzn = im1(tzm);
        
        [rsm, rtm] = eos_s_t(sm, tm, z);
        [rsn, rtn] = eos_s_t(sn, tn, im1(z));
        
        H_imh = ((rsm .* szm - rsn .* szn) ...
            +      (rtm .* tzm - rtn .* tzn)) ...
            ./    ((rsm .* szm + rsn .* szn) ...
            +      (rtm .* tzm + rtn .* tzn) + STRAT_MIN);
        
        % Repeat for 2nd dimension
        [sm, tm] = ppc_val2_mex(Z, SppZ, TppZ, z, 0);
        sn = jm1(sm); tn = jm1(tm);
        [szm, tzm] = ppc_val2_mex(Z, SppZ, TppZ, z, 1);
        szn = jm1(szm); tzn = jm1(tzm);
        
        [rsm, rtm] = eos_s_t(sm, tm, z);
        [rsn, rtn] = eos_s_t(sn, tn, jm1(z));
        
        H_jmh = ((rsm .* szm - rsn .* szn) ...
            +      (rtm .* tzm - rtn .* tzn)) ...
            ./    ((rsm .* szm + rsn .* szn) ...
            +      (rtm .* tzm + rtn .* tzn) + STRAT_MIN);
        
    else
        
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
end


%% Find casts where at least one neighbour provides a valid connection
valid = isfinite(z);
[qu, qt] = bfs_conncomp1(valid, A4, I0);
m = sort(qu(1:qt));

N = length(m);  % number of water columns with at least one valid equation

remap = zeros(ni, nj);  % map from 2D space to 1D vector of water columns
remap(m) = 1 : N;

DXC = OPTS.DXCvec + zeros(ni, nj);
DYC = OPTS.DYCsc  + zeros(ni, nj);

DXCm = DXC(m);
DYCm = DYC(m);
sqrtAREAXm = sqrtAREAX(m);
sqrtAREAYm = sqrtAREAY(m);

%% Matrix
% Build the RHS of the matrix problem
sxm = sx(m);
sym = sy(m);
% 
% good_eqx = ~isnan(sxm) & ~isnan(sxm);
% good_eqy = ~isnan(sym) & ~isnan(sym);
% s2x = sxm(good_eqx) .* sxm(good_eqx);
% s2y = sym(good_eqy) .* sym(good_eqy);

s2xm = sx(m) .* sx(m) .* sqrtAREAXm;
s2ym = sy(m) .* sy(m) .* sqrtAREAYm;

good_eqx = ~isnan(s2xm);
good_eqy = ~isnan(s2ym);
s2x = s2xm(good_eqx);
s2y = s2ym(good_eqy);

rhs = -[s2x; s2y;0];
neqx = length(s2x);
neqy = length(s2y);

neq = neqx + neqy; % Number of equations, excluding density conserving equation.

% Build indices for the rows of the sparse matrix
row = [(1:neq)';(1:neq)'];

% Build indices for the columns of the sparse matrix of slope error
c = remap(A4(:,m).'); 
c_0 = flat(1:N); % z'
c_1 = c(:,1); % jm1
c_2 = c(:,2); % im1

col = [c_2(good_eqx); c_1(good_eqy); c_0(good_eqx); c_0(good_eqy)];

% Build the values of the sparse matrix
p1x = sxm(good_eqx);
p1y = sym(good_eqy);

% val = [-p1x ./ DXCm(good_eqx); -p1y ./ DYCm(good_eqy);...
%         p1x ./ DXCm(good_eqx);  p1y ./ DYCm(good_eqy)];

val = [-p1x ./ DXCm(good_eqx) .* sqrtAREAXm(good_eqx); -p1y ./ DYCm(good_eqy) .* sqrtAREAYm(good_eqy);...
        p1x ./ DXCm(good_eqx) .* sqrtAREAXm(good_eqx);  p1y ./ DYCm(good_eqy) .* sqrtAREAYm(good_eqy)];    
    
val = 2*val;

if STRAT
    % H terms
    H_imh_m = H_imh(m);
    H_jmh_m = H_jmh(m);
    
    vsxjm1 = p1x .* H_imh_m(good_eqx) ./ DXCm(good_eqx) .* sqrtAREAXm(good_eqx);
    vsxj   = p1x .* H_imh_m(good_eqx) ./ DXCm(good_eqx) .* sqrtAREAXm(good_eqx);
    vsyim1 = p1y .* H_jmh_m(good_eqy) ./ DYCm(good_eqy) .* sqrtAREAYm(good_eqy);
    vsyi   = p1y .* H_jmh_m(good_eqy) ./ DYCm(good_eqy) .* sqrtAREAYm(good_eqy);
    
    val2 = [vsxjm1; vsyim1; vsxj; vsyi];
%     val2(~good) = 0;
    val = val + 2*val2;
end

% Build the sparse matrix, with neq+1 rows and nwc columns
mat = sparse(row, col, val, neq, N );

%% New pinning
mr = remap(i0,j0);
% mat(neq+1,mr) = 1e-2;
mat(neq+1,mr) =100* rms(val);

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
if TK > 0
    % using a high-pass filter, specifically a difference operator: grad phi
    tk1x = repelem(1, 1, neqx)';
    tk1y = repelem(1, 1, neqy)';
    val3 = [-tk1x ; -tk1y; tk1x;  tk1y];
    
    % uniform grid of unity
    tik = sparse(row, col, val3, neq, N);
    
    % Apply regularization.  Note tik' * tik is a Laplacian.  Could build that directly...
    mat = mat + TK * (tik' * tik);
    % Note that adding Tikhonov makes pinning not quite perfect.
    % In one test, the solution drifted 6.36cm at ref cast...  current solution is PIN_RIGID
end

%% Solve
sol = mat \ rhs;

dz = zeros(ni, nj);
dz(m) = sol;

end