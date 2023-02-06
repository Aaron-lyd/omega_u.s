function out = rms_Cgrid_v2(T, U, V, AREA_IJ, AREA_iJ, AREA_Ij)
% RMS_CGRID  Calculate Root Mean Square with area weighting appropriate to a C grid
%
% INPUT:
% T: data on the tracer grid
% U: data on the U grid
% V: data on the V grid
% AREA_IJ: area of the tracer grid (often called RAC)
% AREA_iJ: area of the U grid (often called RAW)
% AREA_Ij: area of the V grid (often called RAS).
%
% Soft notation: I, J are integers indexing the tracer grid, and
% i = I - 1/2 and j = J - 1/2 for indexing the U and V grids.
%
% OUTPUT:
% out: the root-mean-square, calculated by summing the non-nan T, U, and V
%      each area weighted by the appropriate area, then dividing by the sum
%      of the areas where the T, U, and V data is non-nan, then finally
%      taking the square root.
%
% EXAMPLE:
% To calculate the rms of |e^sp| + D^f / 1000m, with esp given on the tracer
% grid and slope errors sx and sy on the U and V grids, call:
% rms_Cgrid(esp, sx.^2 * (K / H), sy.^2 * (K / H), AREA_IJ, AREA_iJ, AREA_Ij)
% where H = 1000 % m
%       K = 1000 % m^2 s^-1

AREA_IJ = AREA_IJ + zeros(size(T));
AREA_iJ = AREA_iJ + zeros(size(U));
AREA_Ij = AREA_Ij + zeros(size(V));

flat = @(x) x(:);
sumgood = @(f,g) sum(f(isfinite(g)));  % sum elements of f where g is finite

out = sqrt( ...
  nansum(flat(AREA_IJ .* T.^2)) / sumgood(AREA_IJ, T) + ...
  (nansum(flat(AREA_iJ .* U.^2)) + nansum(flat(AREA_Ij .* V.^2))) ./ (sumgood(AREA_iJ, U) + sumgood(AREA_Ij, V)) );
