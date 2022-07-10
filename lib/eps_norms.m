function [epsL2, epsL1] = eps_norms(s, t, p, use_s_t, wrap, eoses, DIST1_iJ, DIST2_Ij, DIST2_iJ, DIST1_Ij, AREA_iJ, AREA_Ij)
%EPS_NORMS  The L1 and L2 norms of the epsilon neutrality error
%
%
% [epsL2, epsL1] = eps_norms(s, t, p, use_s_t, wrap, eoses, DIST1_iJ, DIST2_Ij, DIST2_iJ, DIST1_Ij, AREA_iJ, AREA_Ij)
% computes the L2 and L1 norms of the neutrality error vector, which itself
% is calculated by ntp_errors.m, on approximately neutral surface on which
% the practical / Absolute salinity is s, the potential / Conservative
% temperature is t, and the pressure or depth is p. If use_s_t is true,
% errors are calculated using s and t differences multiplied by the partial
% derivatives of in-situ density w.r.t s and t, respectively; if use_s_t is
% false, errors are calculated using in-situ density differences and p
% differences multiplied by the partial derivative of in-situ density
% w.r.t. p. Data is treated periodic in the i'th (i=1 or 2) dimension of p
% if and only if wrap(i) is true. The distances and areas of the grid are
% used to weight the neutrality errors when calculating its norm.  See
% Input section for details.
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
%
% For a non-Boussinesq ocean, p is pressure [dbar].  
% For a Boussinesq ocean, p is depth [m, positive].  
%
%
% --- Input:
% s     [ni, nj]: practical / Absolute salinity on the surface
% t     [ni, nj]: potential / Conservative temperature on the surface
% p     [ni, nj]: pressure [dbar] or depth [m, positive] of the surface
% use_s_t [1, 1]: true to compute ex,ey using s and t differences,
%                 false to use eos(s,t,p) and p differences.
% wrap    [2, 1]: wrap(i) is true when the data is periodic in i'th
%                 dimension of p
% DIST1_iJ [ni, nj]: Distance [m] in 1st dimension centred at (I-1/2, J)
% DIST2_Ij [ni, nj]: Distance [m] in 2nd dimension centred at (I, J-1/2)
% DIST2_iJ [ni, nj]: Distance [m] in 2nd dimension centred at (I-1/2, J)
% DIST1_Ij [ni, nj]: Distance [m] in 1st dimension centred at (I, J-1/2)
% AREA_iJ [ni, nj]: Area [m^2] centred at (I-1/2, J). Optional. Should = DIST1_iJ .* DIST2_iJ
% AREA_Ij [ni, nj]: Area [m^2] centred at (I, J-1/2). Optional. Should = DIST1_Ij .* DIST2_Ij
%
% Note: nk is the maximum number of data points per water column,
%       ni is the number of data points in 1st horizontal dimension,
%       nj is the number of data points in 2nd horizontal dimension.
%
% Note: physical units for s,t,p are determined by eos*.m.
%
%
% --- Output:
% epsL2 [ni, nj]: L2 norm of epsilon neutrality error [(kg m^-3) / m or (m^3 / kg) / m]
% epsL1 [ni, nj]: L1 norm of epsilon neutrality error [(kg m^-3) / m or (m^3 / kg) / m]


% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com

flat = @(F) F(:);

% Compute epsilon neutrality errors without handling grid distances
[eps_iJ, eps_Ij] = ntp_errors(s, t, p, 1, 1, use_s_t, false, wrap, eoses);

if nargin < 11
  AREA_iJ = DIST1_iJ .* DIST2_iJ;   % Area [m^2] centred at (I-1/2, J)
  AREA_Ij = DIST1_Ij .* DIST2_Ij;   % Area [m^2] centred at (I, J-1/2)
end

% L2 norm of vector [a_i], weighted by vector [w_i], is sqrt( sum( w_i * a_i^2 ) / sum( w_i ) )
% Here, weights are AREA_iJ and AREA_Ij.
% But also need to divide epsilon by grid distances DIST1_iJ and DIST2_Ij.
% Thus, the numerator of L2 norm needs to multiply epsilon^2 by
%     AREA_iJ ./ DIST1_iJ.^2 = DIST2_iJ ./ DIST1_iJ , 
% and AREA_Ij ./ DIST2_Ij.^2 = DIST1_Ij ./ DIST2_Ij .
epsL2 = sqrt( ...
  (sum(flat( DIST2_iJ ./ DIST1_iJ .* eps_iJ.^2), 'omitnan') + sum(flat( DIST1_Ij ./ DIST2_Ij .* eps_Ij.^2 ), 'omitnan')) ./ ...
  (sum(flat(  AREA_iJ .* isfinite(eps_iJ) ))                + sum(flat(  AREA_Ij .* isfinite(eps_Ij))     ) ) ...
  );

if nargout < 2; return; end

% L1 norm of vector [a_i], weighted by vector [w_i], is sum( w_i * |a_i| ) / sum( w_i )
% Here, weights are AREA_iJ and AREA_Ij.
% But also need to divide epsilon by grid distances DIST1_iJ and DIST2_Ij.
% Thus, the numerator of L1 norm needs to multiply epsilon by
%     AREA_iJ ./ DIST1_iJ = DIST2_iJ , 
% and AREA_Ij ./ DIST2_Ij = DIST1_Ij .
epsL1 =  ...
  (sum(flat( DIST2_iJ .* abs(eps_iJ) ), 'omitnan' ) + sum(flat( DIST1_Ij .* abs(eps_Ij) ), 'omitnan' )) ./ ...
  (sum(flat(  AREA_iJ .* isfinite(eps_iJ) ))        + sum(flat(  AREA_Ij .* isfinite(eps_Ij) ) ) ) ;
