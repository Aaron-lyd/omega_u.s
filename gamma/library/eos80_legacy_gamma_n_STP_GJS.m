function [gamma_n, gamma_error_lower, gamma_error_upper, SP_ntp, t_ntp, T_ntp, p_ntp,...
          gamma_ntp1, gamma_ntp2, gamma_ntp3, gamma_ntp4,...
          SP_ntp1, t_ntp1, T_ntp1, p_ntp1, SP_ntp2, t_ntp2, T_ntp2, p_ntp2,...
          SP_ntp3, t_ntp3, T_ntp3, p_ntp3, SP_ntp4, t_ntp4, T_ntp4, p_ntp4] = eos80_legacy_gamma_n_STP_GJS(SP,t,p,long,lat)
% Also output salinity, temperature, and pressure of the bottle where gamma_n is taken
% Modification by Geoff Stanley
% 13/05/2019.
% Bug fix 7/6/2019: Eliminated usage of "if <not a scalar quantity>"

% eos80_legacy_gamma_n                                      neutral density
%==========================================================================
%
% USAGE:  
%  [gamma_n,{gamma_error_lower,gamma_error_upper}] = ...
%                             eos80_legacy_gamma_n(SP,t,p,long,lat)
%
% DESCRIPTION:
%  Calculates neutral density value, gamma_n, in the open ocean by 
%  spatially interpolating the global reference data set, gamma_n_ref, to
%  the location of the seawater sample.  The input variables are those of
%  EOS-80, namely Practical Salinity and in-situ temperature.
% 
%  This function uses the gamma_n look up table of Jackett & McDougall 
%  (1997). 
%
% INPUT:
%  SP    =  Practical Salinity                                 [ unitless ]
%  t     =  in-situ temperature (ITS-90)                          [ deg C ]
%  p     =  sea pressure                                           [ dbar ] 
%          ( i.e. absolute pressure - 10.1325 dbar )
%  long  =  longitude in decimal degrees                     [ 0 ... +360 ]
%                                                      or [ -180 ... +180 ]
%  lat   =  latitude in decimal degrees north               [ -90 ... +90 ]
%
%  SP and t need to have the same dimensions.
%  p, lat & long may have dimensions 1x1 or Mx1 or 1xN or MxN,
%  where SP & t is MxN.
%
% OUTPUT:
%  gamma_n  =  neutral density                                  [ kg m^-3 ]
% Optional:
%  gamma_error_lower = lower estimate of the error
%  gamma_error_upper = upper estimate of the error
%
% AUTHOR: 
%  David Jackett                                       [ help@teos-10.org ]
%
% MODIFIED:
%  Paul Barker and Trevor McDougall
%
% VERSION NUMBER: 3.05.10 (6th March, 2018)
%
% REFERENCES:
% Jackett, D.R., and T.J. McDougall, 1997: A Neutral Density Variable for
%  the World's Oceans. J. Phys. Oceanogr., 27, 237�263. 
%  doi: 10.1175/1520-0485(1997)0272.0.CO;2
%
% Unesco, 1983: Algorithms for computation of fundamental properties of 
%  seawater. Unesco technical papers in marine science 44, 53pp.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Check variables and resize if necessary
%--------------------------------------------------------------------------

if ~(nargin == 5)
   error('eos80_legacy_gamma_n:  Requires five inputs')
end

sz_SP = size(SP);
sz_t = size(t);

if sum(sz_SP - sz_t) ~= 0 
    error('eos80_legacy_gamma_n: SP and t need be of the same size')
end 

p = gsw_resize(p,sz_SP);
long = gsw_resize(long,sz_SP);
lat = gsw_resize(lat,sz_SP);

long(long < 0) = long(long < 0) + 360; 

if any(p < -1.5)
    error('eos80_legacy_gamma_n: pressure needs to be positive')
end

%set any pressures between 0 and 1.5 to be equal to 0 (i.e. the surface)
p(p < 0) = 0;

%{
if any(long < 0 | long > 360 | lat < -80 | lat > 64)
    [Ioor] = find(long < 0 | long > 360 | lat < -80 | lat > 64);
    SP(Ioor) = NaN;
end

if any(SP(:) < 0 | SP(:) > 42 | t(:) < -2.5 | t(:) > 40 | p(:) < 0 | p(:) > 10000)
    [Ioor] = find(SP(:) < 0 | SP(:) > 42 | t(:) < -2.5 | t(:) > 40 | p(:) < 0 | p(:) > 10000);
    SP(Ioor) = NaN;
end
%}
% GJS bug fix 07/06/2019.   The above original code only tests the first
% column output by the any function.
SP(long < 0 | long > 360 | lat < -80 | lat > 64 | SP < 0 | SP > 42 | t < -2.5 | t > 40 | p < 0 | p > 10000) = NaN;


Iocean = ~isnan(SP + t + p + lat + long);

gamma_n = nan(size(SP));
SP_ntp = nan(size(SP));
t_ntp = nan(size(SP));
T_ntp = nan(size(SP));
p_ntp = nan(size(SP));

SP_ntp1 = nan(size(SP));
SP_ntp2 = nan(size(SP));
SP_ntp3 = nan(size(SP));
SP_ntp4 = nan(size(SP));

t_ntp1 = nan(size(SP));
t_ntp2 = nan(size(SP));
t_ntp3 = nan(size(SP));
t_ntp4 = nan(size(SP));

T_ntp1 = nan(size(SP));
T_ntp2 = nan(size(SP));
T_ntp3 = nan(size(SP));
T_ntp4 = nan(size(SP));

p_ntp1 = nan(size(SP));
p_ntp2 = nan(size(SP));
p_ntp3 = nan(size(SP));
p_ntp4 = nan(size(SP));

gamma_ntp1 = nan(size(SP));
gamma_ntp2 = nan(size(SP));
gamma_ntp3 = nan(size(SP));
gamma_ntp4 = nan(size(SP));



if nargout > 1
    gamma_error_lower = gamma_n;
    gamma_error_upper = gamma_n;
    [gamma_n(Iocean), gamma_error_lower(Iocean), gamma_error_upper(Iocean), SP_ntp(Iocean), t_ntp(Iocean), T_ntp(Iocean), p_ntp(Iocean),...
          gamma_ntp1(Iocean), gamma_ntp2(Iocean), gamma_ntp3(Iocean), gamma_ntp4(Iocean),...
          SP_ntp1(Iocean), t_ntp1(Iocean), T_ntp1(Iocean), p_ntp1(Iocean), SP_ntp2(Iocean), t_ntp2(Iocean), T_ntp2(Iocean), p_ntp2(Iocean),...
          SP_ntp3(Iocean), t_ntp3(Iocean), T_ntp3(Iocean), p_ntp3(Iocean), SP_ntp4(Iocean), t_ntp4(Iocean), T_ntp4(Iocean), p_ntp4(Iocean)] = ...
        eos80_legacy_gamma_n_labelling_STP_GJS(SP(Iocean),t(Iocean),p(Iocean),long(Iocean),lat(Iocean));
else
    [gamma_n(Iocean), ~, ~, SP_ntp(Iocean), t_ntp(Iocean), T_ntp(Iocean), p_ntp(Iocean),...
          gamma_ntp1(Iocean), gamma_ntp2(Iocean), gamma_ntp3(Iocean), gamma_ntp4(Iocean),...
          SP_ntp1(Iocean), t_ntp1(Iocean), T_ntp1(Iocean), p_ntp1(Iocean), SP_ntp2(Iocean), t_ntp2(Iocean), T_ntp2(Iocean), p_ntp2(Iocean),...
          SP_ntp3(Iocean), t_ntp3(Iocean), T_ntp3(Iocean), p_ntp3(Iocean), SP_ntp4(Iocean), t_ntp4(Iocean), T_ntp4(Iocean), p_ntp4(Iocean)] = ...
        eos80_legacy_gamma_n_labelling_STP_GJS(SP(Iocean),t(Iocean),p(Iocean),long(Iocean),lat(Iocean));
end


end
