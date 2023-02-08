function [SP_ntp1, t_ntp1, p_ntp1, theta_ntp1,...
          SP_ntp2, t_ntp2, p_ntp2, theta_ntp2,...
          SP_ntp3, t_ntp3, p_ntp3, theta_ntp3,...
          SP_ntp4, t_ntp4, p_ntp4, theta_ntp4] = eos80_gamma_scv_stp(SP,t,p,long,lat)

% eos80_legacy_gamma_n_labelling                            neutral density
%==========================================================================
%
% USAGE:
%  [gamma_n, gamma_error_lower, gamma_error_upper] = ...
%                    eos80_legacy_gamma_n_labelling(SP,t,p,long,lat)
%
% DESCRIPTION:
%  Calculates neutral density value, gamma_n, in the open ocean by 
%  spatially interpolating the global reference data set, gamma_n_ref, to
%  the location of the seawater sample using the neutral tangent plane.
%  
%  The correct way of calculating gamma_n is by calling 
%  eos80_legacy_gamma_n.  
%
%  This function uses version 1.0 of the gamma_n look up table (1997). 
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
%  SP, t, p, long & lat need to be vectors and have the same dimensions.
%
% OUTPUT:
%  gamma_n  =  neutral density                                  [ kg m^-3 ]
%  SP_ntp = the practical salinity on the reference cast       [ unitless ]
%  T_ntp = the in-situ temperature on the reference cast          [ deg C ]
%  P_ntp = the sea pressure on the reference cast                  [ dbar ]
% Optional:
%  gamma_error_lower = lower estimate of the error
%  gamma_error_upper = upper estimate of the error
%
% AUTHOR: 
%  David Jackett                                       [ help@teos-10.org ]
%
% MODIFIED:
%  Yandong Lang (08/04/2019)
%
% VERSION NUMBER: 3.05.10 (6th March, 2018)
%
% REFERENCES:
%  Jackett, D.R., and T.J. McDougall, 1997: A Neutral Density Variable for
%   the World's Oceans. J. Phys. Oceanogr., 27, pp. 237–263.
%
%  The software is available from http://www.TEOS-10.org
%
%==========================================================================

%--------------------------------------------------------------------------
% Start of the calculation (extracting from a look up table)
%--------------------------------------------------------------------------
persistent gamma_n_ref lat_ref long_ref p_ref SP_ref t_ref a_ref

if isempty(gamma_n_ref)
    gamma_n_data = 'eos80_legacy_gamma_n_data.mat';
    
    gamma_n_data_file = which(gamma_n_data);

    load (gamma_n_data_file,'gamma_n_ref','lat_ref','long_ref','p_ref',...
        'a_ref','SP_ref','t_ref','n_ref');
end

% precalculate constants 
nx = length(long_ref); 
ny = length(lat_ref); 
nz = length(p_ref); 

% Calculate the 4 grid points surrounding the data
Ix0 = floor(1 + (nx-1)*(long - long_ref(1))./(long_ref(nx) - long_ref(1)));
Ix0 = Ix0(:); 
Ix0(Ix0 == nx) = nx - 1;
              
Iy0 = floor(1 + (ny-1)*(lat - lat_ref(1))./(lat_ref(ny) - lat_ref(1)));
Iy0 = Iy0(:); 
Iy0(Iy0 == ny) = ny - 1;
     
Iy0_Ix0ny = Iy0 + Ix0.*ny;        
I1 = Iy0_Ix0ny - ny;              
I2 = Iy0_Ix0ny;
I3 = Iy0_Ix0ny + 1;
I4 = Iy0_Ix0ny + (1 - ny);

% calculate the ntp for the bottle to each of the 4 surrounding grid points
[mp,np] = size(I1);
p_ref1_dummy = repmat(p_ref(:),1,mp);
p_ref1_dummy(isnan(SP_ref(:,I1))) = NaN;
[mp,np] = size(I2);
p_ref2_dummy = repmat(p_ref(:),1,mp);
p_ref2_dummy(isnan(SP_ref(:,I2))) = NaN;
[mp,np] = size(I3);
p_ref3_dummy = repmat(p_ref(:),1,mp);
p_ref3_dummy(isnan(SP_ref(:,I3))) = NaN;
[mp,np] = size(I4);
p_ref4_dummy = repmat(p_ref(:),1,mp);
p_ref4_dummy(isnan(SP_ref(:,I4))) = NaN;

% do the quick vectorised ntp calculation
[SP_ntp1, t_ntp1, p_ntp1] = eos80_legacy_scv_bottle_to_cast(SP(:).',t(:).',p(:).',SP_ref(:,I1),t_ref(:,I1),p_ref1_dummy);
[SP_ntp2, t_ntp2, p_ntp2] = eos80_legacy_scv_bottle_to_cast(SP(:).',t(:).',p(:).',SP_ref(:,I2),t_ref(:,I2),p_ref2_dummy);
[SP_ntp3, t_ntp3, p_ntp3] = eos80_legacy_scv_bottle_to_cast(SP(:).',t(:).',p(:).',SP_ref(:,I3),t_ref(:,I3),p_ref3_dummy);
[SP_ntp4, t_ntp4, p_ntp4] = eos80_legacy_scv_bottle_to_cast(SP(:).',t(:).',p(:).',SP_ref(:,I4),t_ref(:,I4),p_ref4_dummy);

% % calculate spatial weights, based on their distance from each grid point
% rx = (long(:).' - long_ref(Ix0))./(long_ref(Ix0+1) - long_ref(Ix0));
% ry = (lat(:).' - lat_ref(Iy0))./(lat_ref(Iy0+1) - lat_ref(Iy0));
% 
% SP_ntp = (1-ry).*(SP_ntp1 + rx.*(SP_ntp2 - SP_ntp1)) + ry.*(SP_ntp4 + rx.*(SP_ntp3 - SP_ntp4));
% t_ntp = (1-ry).*(t_ntp1 + rx.*(t_ntp2 - t_ntp1)) + ry.*(t_ntp4 + rx.*(t_ntp3 - t_ntp4));
% p_ntp = (1-ry).*(p_ntp1 + rx.*(p_ntp2 - p_ntp1)) + ry.*(p_ntp4 + rx.*(p_ntp3 - p_ntp4));

% extrap_weights = ones(4,length(gamma_ntp1));
% % replace any NaN's with the mean of the non-NaN grid points
% if any(isnan(sum(gamma_ntps))) | any(sum(extrap_weights) < 4)
%     [I_nan] = find(isnan(sum(gamma_ntps)) | sum(extrap_weights) < 1);
%     gamma_ntps(:,I_nan) = gsw_refdata_replace_nan(gamma_ntps(:,I_nan), extrap_weights(:,I_nan));
% end

% SP, t, p on the reference cast
SP_ntp1 = reshape(SP_ntp1, size(SP));
SP_ntp2 = reshape(SP_ntp2, size(SP));
SP_ntp3 = reshape(SP_ntp3, size(SP));
SP_ntp4 = reshape(SP_ntp4, size(SP));
% SP_ntp = reshape(SP_ntp, size(SP));

SP_ntp1 = double(SP_ntp1);
SP_ntp2 = double(SP_ntp2);
SP_ntp3 = double(SP_ntp3);
SP_ntp4 = double(SP_ntp4);
% SP_ntp = double(SP_ntp);

t_ntp1 = reshape(t_ntp1, size(SP));
t_ntp2 = reshape(t_ntp2, size(SP));
t_ntp3 = reshape(t_ntp3, size(SP));
t_ntp4 = reshape(t_ntp4, size(SP));
% t_ntp = reshape(t_ntp, size(SP));

t_ntp1 = double(t_ntp1);
t_ntp2 = double(t_ntp2);
t_ntp3 = double(t_ntp3);
t_ntp4 = double(t_ntp4);
% t_ntp = double(t_ntp);

p_ntp1 = reshape(p_ntp1, size(SP));
p_ntp2 = reshape(p_ntp2, size(SP));
p_ntp3 = reshape(p_ntp3, size(SP));
p_ntp4 = reshape(p_ntp4, size(SP));
% p_ntp = reshape(p_ntp, size(SP));

p_ntp1 = double(p_ntp1);
p_ntp2 = double(p_ntp2);
p_ntp3 = double(p_ntp3);
p_ntp4 = double(p_ntp4);
% p_ntp = double(p_ntp);


% SA_ntp1 = gsw_SA_from_SP(SP_ntp1, p_ntp1, long, lat);
% SA_ntp2 = gsw_SA_from_SP(SP_ntp2, p_ntp2, long, lat);
% SA_ntp3 = gsw_SA_from_SP(SP_ntp3, p_ntp3, long, lat);
% SA_ntp4 = gsw_SA_from_SP(SP_ntp4, p_ntp4, long, lat);
% 
% SA_ntp1 = double(SA_ntp1);
% SA_ntp2 = double(SA_ntp2);
% SA_ntp3 = double(SA_ntp3);
% SA_ntp4 = double(SA_ntp4);

theta_ntp1 = eos80_legacy_pt(SP_ntp1,t_ntp1,p_ntp1,zeros(size(SP)));
theta_ntp2 = eos80_legacy_pt(SP_ntp2,t_ntp2,p_ntp2,zeros(size(SP)));
theta_ntp3 = eos80_legacy_pt(SP_ntp3,t_ntp3,p_ntp3,zeros(size(SP)));
theta_ntp4 = eos80_legacy_pt(SP_ntp4,t_ntp4,p_ntp4,zeros(size(SP)));

theta_ntp1 = reshape(theta_ntp1, size(SP));
theta_ntp2 = reshape(theta_ntp2, size(SP));
theta_ntp3 = reshape(theta_ntp3, size(SP));
theta_ntp4 = reshape(theta_ntp4, size(SP));

% theta_ntps(1,:) = theta_ntp1(:);
% theta_ntps(2,:) = theta_ntp2(:);
% theta_ntps(3,:) = theta_ntp3(:);
% theta_ntps(4,:) = theta_ntp4(:);
% 
% theta_ntp = (1-ry).*(theta_ntps(1,:) + rx.*(theta_ntps(2,:) - theta_ntps(1,:))) + ry.*(theta_ntps(4,:) + rx.*(theta_ntps(3,:) - theta_ntps(4,:)));
% theta_ntp = reshape(theta_ntp, size(SP));

theta_ntp1 = double(theta_ntp1);
theta_ntp2 = double(theta_ntp2);
theta_ntp3 = double(theta_ntp3);
theta_ntp4 = double(theta_ntp4);
% theta_ntp = double(theta_ntp);


end