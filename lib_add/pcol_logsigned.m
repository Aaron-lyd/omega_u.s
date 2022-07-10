function [hp, hc] = pcol_logsigned(ax, x, y, z, m, M, s, col_m, col_n, cbarpos)
% PCOL_LOGSIGNED  A signed version of mapping log(abs(z)).
%
% INPUT:
%   ax      axes handle to plot into
%   x       [], 1D, or 2D: x-axis
%   y       [], 1D, or 2D: y-axis
%   z       variable whose signed log gets plotted in color
%   m       integer: |z| restricted to minimum value of 10^m
%   M       integer: |z| restricted to maximum value of 10^M
%   s       integer: number of colors between consecutive powers of 10
% col_m     3 element array: RGB values for values with |z| < 10^m
% col_n     3 element array: RGB values for nan values
% cbarpos   If 4 element array: left, bottom, width, height values for color bar position.
%           If false, no colorbar is added. 
%           If omitted or [], vertical colorbar is added in the default place. 
%
% OUTPUT:
%   hp      handle to pcolor
%   hc      handle to colorbar
%
% EXAMPLE:
% dat = [-150, -50, -.2, 0, .8, 1.01, 12, -5,  98, 102]; 
% figure; ax=axes; pcol_logsigned(ax, 1, 1:length(dat), dat, 0, 2);
%
% Requirements: cbrewer, pcol
%
% Author: Geoff Stanley.
% Acknowledgements:  inspired by plot_grid_logscale.m by Graeme MacGilchrist

if nargin < 7 || isempty(s)
  s = 1;
end

n = 2 * (M - m) * s + 1; % # of colors:  s * (# of intervals, positive and negative), plus 1 for middle color. 
cm = flipud(cbrewer2('div', 'RdBu', n)); % (negative) Blue -> White -> Red (positive)

if nargin < 8 || isempty(col_m)
  col_m = cm((n+1)/2,:); % pale white
end

if nargin < 9 || isempty(col_n)
  col_n = [1 1 1] * .75; % light grey
end

% Repeat zero color in the middle, s times.  After transforming z (next),
% values between -0.5 and 0.5 will be mapped to the zero color.  (The only
% value z will take between -0.5 and 0.5 will be, identically, 0.)
cm = [ ...
  cm(1 : (n-1)/2,:); ...   % blue -> white
  repmat(col_m, s, 1); ... % middle colour
  cm((n+3)/2 : n,:)];      % white -> red

% First, restrict log10(|z|) to the range [m, M], then shift to the range
% [0, M-m] and reimpose the sign of z.  Any elements with log10(|z|) < m
% will now be 0.
% Second, shift to the range [1/2, (M - m) + 1/2], and reimpose the sign of
% z.
% Note that sign(0) == 0, so again elements with log10(|z|) < m will remain
% 0.
% Final result: z is in [-(M - m) - 1/2, -1/2] U {0} U [1/2, (M - m) + 1/2]
z = sign(z) .* max(min(log10(abs(z)) - m, M-m), 0);
z = sign(z) .* (abs(z) + 0.5);

% Make pcolor plot
hp = pcol(ax, x, y, z, col_n);
ax.CLim = [-(M-m)-.5, M-m+.5];
colormap(ax, cm);

% Add colorbar
if nargin == 10 && islogical(cbarpos) && cbarpos == false
  return
end

if nargin == 10 && length(cbarpos) == 4
  w = cbarpos(3);
  h = cbarpos(4);
  if h > w
    % vertical colorbar (location east, then override)
    hc = colorbar('Location', 'eastoutside', 'Position', cbarpos);
  else
    % horizontal colorbar (location south, then override)
    hc = colorbar('Location', 'southoutside', 'Position', cbarpos);
  end
else
  hc = colorbar;
end

hc.YTick = ax.CLim(1) : ax.CLim(2);

hc.YTickLabel = [ ...
  arrayfun(@(F) ['-10^{' num2str(F) '}'], M : -1 : m, 'UniformOutput', false), ...
  arrayfun(@(F) [ '10^{' num2str(F) '}'], m :  1 : M, 'UniformOutput', false) ] ;

hc.TickLength = 0.015; % longer ticks
