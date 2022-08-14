function [s, s2, s_rms, s2_rms, e_hel, e_hel_rms, flux, d] = s_ehel_flux(s, t, z, S, T, Z, U, V, wrap, OPT)

% functions
im1 = @(F) circshift(F, [+1, 0]); % if x is [lat x lon], this shifts south
jm1 = @(F) circshift(F, [0, +1]); % if x is [lat x lon], this shifts west

ip1 = @(D) circshift(D, [-1 0]);
jp1 = @(D) circshift(D, [0 -1]);

A_X_f = @(F) (F + ip1(F)) / 2;
A_Y_f = @(F) (F + jp1(F)) / 2;

lead1 = @(x) reshape(x, [1 size(x)]);
flat = @(F) F(:);

[nx, ny] = size(s);

%% Grid information
weight_u = OPT.RAWvec + zeros(OPT.nx, OPT.ny);
weight_v = OPT.RASvec + zeros(OPT.nx, OPT.ny);
dx = OPT.DXCvec;
dy = OPT.DYCsc;

%% Interpolation
% interpfn = @ppc_linterp;
interpfn = @ppc_pchip;
SppZ = interpfn(Z, S);
TppZ = interpfn(Z, T);
UppZ = interpfn(Z, U); % knowing Z is 1D vector, horizontally uniform
VppZ = interpfn(Z, V); % knowing Z is 1D vector, horizontally uniform

DATACUBE = OPT.data_cube;

%% slope error
mode_s = OPT.MODES;
if mode_s == 1 % for omege+
    [~,~,sx,sy] = ntp_errors(s, t, z, dx, dy, 1, 0, wrap,{}, 9.81, S, T, Z);
else
    tolz = OPT.tolz;
    [sx, sy] = ntp_slope_error(SppZ, TppZ, Z, z, tolz, dx, dy, wrap);
end

% % Method 1 to calculate s2 and s
% s2 = (sx.^2 + (circshift(sx, [-1,0])).^2 + sy.^2 + (circshift(sy, [0,-1])).^2)/2;
% s2_rms = s2 - nanmean(s2(:));
% good = ~isnan(s2_rms);
% s2_rms = rms(s2_rms(good));
% 
% s = sqrt(s2);
% s_rms = s - nanmean(s(:));
% good = ~isnan(s_rms);
% s_rms = rms(s_rms(good));

% Method 2 to calculate s2 and s
sx2 = sx .* sx .* (OPT.DXCvec .* OPT.DYGsc);
sy2 = sy .* sy .* (OPT.DXGvec .* OPT.DYCsc);
s2 = (sx2 + ip1(sx2) + sy2 + jp1(sy2)) ./ (2 * OPT.RACvec);
s = sqrt(s2);

% s2_rms = s2 - nanmean(s2(:));
% good = ~isnan(s2_rms);
% s2_rms = rms(s2_rms(good));
% 
% s_rms = s - nanmean(s(:));
% good = ~isnan(s_rms);
% s_rms = rms(s_rms(good));
% 
% sx_rms = sx - nanmean(sx(:));
% good = ~isnan(sx_rms);
% sx_rms = rms(sx_rms(good));
% 
% sy_rms = sy - nanmean(sy(:));
% good = ~isnan(sy_rms);
% sy_rms = rms(sy_rms(good));

% s2_rms = s2 - nanmean(s2(:));
% good = ~isnan(s2_rms);
% s2_rms = sqrt(sum(flat(s2_rms(good).^2))/sum(good(:)));
% 
% s_rms = s - nanmean(s(:));
% good = ~isnan(s_rms);
% s_rms = sqrt(sum(flat(s_rms(good).^2))/sum(good(:)));
% 
% sx_rms = sx - nanmean(sx(:));
% good = ~isnan(sx_rms);
% sx_rms = rms(sx_rms(good));
% 
% sy_rms = sy - nanmean(sy(:));
% good = ~isnan(sy_rms);
% sy_rms = rms(sy_rms(good));

goodx = ~isnan(sx);
goody = ~isnan(sy);

if isscalar(OPT.DXCvec)
    AreaX = OPT.DXCvec .* OPT.DYGsc;
    AreaY = OPT.DXGvec .* OPT.DYCsc;
    s2x_grid = sx(goodx) .^2 .* AreaX;
    s2y_grid = sy(goody) .^2 .* AreaY;
    s4x_grid = sx(goodx) .^4 .* AreaX;
    s4y_grid = sy(goody) .^4 .* AreaY;
else
    
    AreaX = repmat((OPT.DXCvec .* OPT.DYGsc), [nx 1]);
    AreaY = repmat((OPT.DXGvec .* OPT.DYCsc), [nx 1]);
    
    s2x_grid = sx(goodx) .^2 .* AreaX(goodx);
    s2y_grid = sy(goody) .^2 .* AreaY(goody);
    s4x_grid = sx(goodx) .^4 .* AreaX(goodx);
    s4y_grid = sy(goody) .^4 .* AreaY(goody);
end

sx_rms = sqrt(nansum(s2x_grid(:))/nansum(AreaX(:)));
sy_rms = sqrt(nansum(s2y_grid(:))/nansum(AreaY(:)));

s_rms = sqrt( (nansum(s2x_grid(:))+nansum(s2y_grid(:)))/(nansum(AreaX(:))+nansum(AreaY(:))));
s2_rms = sqrt( (nansum(s4x_grid(:))+nansum(s4y_grid(:)))/(nansum(AreaX(:))+nansum(AreaY(:))));


if ~DATACUBE
    sx_st = sx(:,136:146);
    sy_st = sy(:,136:146);
end

goodx_st = ~isnan(sx_st);
goody_st = ~isnan(sy_st);


if isscalar(OPT.DXCvec)
    AreaX = OPT.DXCvec .* OPT.DYGsc;
    AreaY = OPT.DXGvec .* OPT.DYCsc;
    s2stx_grid = sx_st(goodx_st) .^2 .* AreaX;
    s2sty_grid = sy_st(goody_st) .^2 .* AreaY;
    s4stx_grid = sx_st(goodx_st) .^4 .* AreaX;
    s4sty_grid = sy_st(goody_st) .^4 .* AreaY;
else
    
    AreaX = repmat((OPT.DXCvec .* OPT.DYGsc), [nx 1]);
    AreaY = repmat((OPT.DXGvec .* OPT.DYCsc), [nx 1]);
    if ~DATACUBE
        AreaX_st = AreaX(:,136:146);
        AreaY_st = AreaY(:,136:146);
        
        s2stx_grid = sx_st(goodx_st) .^2 .* AreaX_st(goodx_st);
        s2sty_grid = sy_st(goody_st) .^2 .* AreaY_st(goody_st);
        s4stx_grid = sx_st(goodx_st) .^4 .* AreaX_st(goodx_st);
        s4sty_grid = sy_st(goody_st) .^4 .* AreaY_st(goody_st);
    end
end

if ~DATACUBE
    sx_st_rms = sqrt(nansum(s2stx_grid(:))/nansum(AreaX_st(:)));
    sy_st_rms = sqrt(nansum(s2sty_grid(:))/nansum(AreaY_st(:)));
    
    s_st_rms = sqrt( (nansum(s2stx_grid(:))+nansum(s2sty_grid(:)))/(nansum(AreaX_st(:))+nansum(AreaY_st(:))));
    s2_st_rms = sqrt( (nansum(s4stx_grid(:))+nansum(s4sty_grid(:)))/(nansum(AreaX_st(:))+nansum(AreaY_st(:))));
end

%% ehel and flux
% u, v on the u and v grid
u = ppc_val(Z, UppZ, (z + im1(z)) / 2, 0);
v = ppc_val(Z, VppZ, (z + jm1(z)) / 2, 0);

u2 = u .* u .* (OPT.DXCvec .* OPT.DYGsc);
v2 = v .* v .* (OPT.DXGvec .* OPT.DYCsc);
u_t2 = (u2 + ip1(u2) + v2 + jp1(v2)) ./ (2 * OPT.RACvec);
u_t = sqrt(u_t2); % velocity on the tracer cell

% Calculate e^hel on (i-1/2, j):
e_hel_u = u.*sx;
% Calculate e^hel on (i, j-1/2):
e_hel_v = v.*sy;

e_hel_u_rms = e_hel_u - nanmean(e_hel_u(:));
good = ~isnan(e_hel_u_rms);
e_hel_u_rms = rms(e_hel_u_rms(good));

e_hel_v_rms = e_hel_v - nanmean(e_hel_v(:));
good = ~isnan(e_hel_v_rms);
e_hel_v_rms = rms(e_hel_v_rms(good));

mode_hel = OPT.MODEHEL;
if mode_hel == 1 % divergence method to calculate ehel
    
    ehelx = u .* sx .* (OPT.DXCvec .* OPT.DYGsc);
    ehely = v .* sy .* (OPT.DXGvec .* OPT.DYCsc);
    e_hel = (ehelx + ip1(ehelx) + ehely + jp1(ehely)) ./ (2 * OPT.RACvec);
    
    flux = e_hel .* OPT.RACvec;
    flux = nansum(flux(:));
    
    e_hel_rms = e_hel - nanmean(e_hel(:));
    good = ~isnan(e_hel_rms);
    e_hel_rms = rms(e_hel_rms(good));
    
    if ~DATACUBE
        
        e_hel_st = e_hel(:,136:146);
        
        e_hel_st_rms = e_hel_st - nanmean(e_hel_st(:));
        good = ~isnan(e_hel_st_rms);
        e_hel_st_rms = rms(e_hel_st_rms(good));
    end
    
    flux_u = e_hel_u .* weight_u;
    flux_v = e_hel_v .* weight_v;
    
    flux_u = nansum(flux_u(:));
    flux_v = nansum(flux_v(:));
else
    e_hel = e_hel_u + e_hel_v;
    
    % e_hel_rms
    e_hel_rms = e_hel - nanmean(e_hel(:));
    good = ~isnan(e_hel_rms);
    e_hel_rms = rms(e_hel_rms(good));
    
    if ~DATACUBE
        
        e_hel_st = e_hel(:,136:146);
        
        e_hel_st_rms = e_hel_st - nanmean(e_hel_st(:));
        good = ~isnan(e_hel_st_rms);
        e_hel_st_rms = rms(e_hel_st_rms(good));
    end
    
    % flux
    flux_u = e_hel_u .* weight_u;
    flux_v = e_hel_v .* weight_v;
    
    flux_u = nansum(flux_u(:));
    flux_v = nansum(flux_v(:));
    
    flux = flux_u + flux_v;


    flux = flux_u + flux_v;
end

flux_u_positive = nansum((e_hel_u(e_hel_u(:)>0)).*(weight_u(e_hel_u(:)>0)));
flux_u_negative = nansum((e_hel_u(e_hel_u(:)<0)).*(weight_u(e_hel_u(:)<0)));
flux_v_positive = nansum((e_hel_v(e_hel_v(:)>0)).*(weight_v(e_hel_v(:)>0)));
flux_v_negative = nansum((e_hel_v(e_hel_v(:)<0)).*(weight_v(e_hel_v(:)<0)));

%% detail output for diagnostics
d.u = u;
d.v = v;
d.u_meanT = A_X_f(u);
d.v_meanT = A_Y_f(v);
d.u_t = u_t;
d.sx = sx;
d.sy = sy;
d.sx_rms = sx_rms;
d.sy_rms = sy_rms;
if ~DATACUBE
    d.s2_st_rms = s2_st_rms;
end
d.e_hel_u = e_hel_u;
d.e_hel_v = e_hel_v;
d.e_hel_u_rms = e_hel_u_rms;
d.e_hel_v_rms = e_hel_v_rms;
if ~DATACUBE
    d.e_hel_st_rms = e_hel_st_rms;
end
d.flux_u = flux_u;
d.flux_v = flux_v;
d.flux_u_positive = flux_u_positive;
d.flux_v_positive = flux_v_positive;
d.flux_u_negative = flux_u_negative;
d.flux_v_negative = flux_v_negative;

end



