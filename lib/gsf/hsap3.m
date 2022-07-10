function Y = hsap3(P, ATMP, ETAN, A, grav, varargin)
%HSAP3  Integrate hydrostatic balance to obtain acceleration potential at every data point.
%
%
% Y = hsap3(P, ATMP, ETAN, A, grav)
% integrates hydrostatic balance to obtain the acceleration potential Y --
% the depth times the gravitational acceleration -- at the data sites of A
% and P, in an ocean with atmospheric pressure ATMP, sea-surface height
% ETAN, specific volume A at the data sites of pressure P, and
% gravitational acceleration grav.
%
% Y = hsap3(Z, ATMP, ETAN, R, grav, rho_c)
% integrates hydrostatic balance to obtain the acceleration potential Y --
% the pressure divided by the Boussinesq reference density -- at the data
% sites of R and Z, in a Boussinesq ocean with atmospheric pressure ATMP,
% sea-surface height ETAN, in-situ density R at the data sites of depth Z,
% with gravitational acceleration grav, and Boussinesq reference density
% rho_c.
%
%
% --- Input:
% ATMP [ni, nj]: atmospheric pressure loading [dbar]
% ETAN [ni, nj]: sea-surface height [m]
% A [nk, ni, nj]: specific volume [m^3 kg^-1]
% R [nk, ni, nj]: in-situ density [kg m^-3]
% Z [nk, ni, nj] or [nk, 1]: depth [m, positive]
% P [nk, ni, nj] or [nk, 1]: pressure [dbar]
% grav [1, 1]: gravitational acceleration [m s^-2]
% rho_c [1, 1]: Boussinesq reference density [kg m^-3]
%
%
% --- Output:
% Y [nk, ni, nj]: acceleration potential from hydrostatic balance [m^2 s^-2]

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com



narginchk(5,6);

db2Pa = 1e4;
lead1 = @(p) reshape(p, [1 size(p)]);

if nargin == 6
    % Boussinesq form
    rho_c = varargin{1};
    if isvector(P)
        Y = (db2Pa/rho_c) * lead1(ATMP) + grav * lead1(ETAN) ...
            + cumsum((grav/(2*rho_c)) * diff([0; P]) .* (A([1, 1:end-1],:,:) + A), 1);
    else
        Y = (db2Pa/rho_c) * lead1(ATMP) + grav * lead1(ETAN) ...
            + cumsum((grav/(2*rho_c)) * cat(1, P(1,:,:), diff(P,[],1)) .* (A([1, 1:end-1],:,:) + A), 1);
    end
    
else
    % Non-Boussinesq form
    Y = + grav * lead1(ETAN) ...
        - (db2Pa / 2) * cumsum(cat(1, P(1,:,:) - lead1(ATMP), diff(P,[],1)) .* (A([1 1:end-1],:,:) + A), 1);
    
end
