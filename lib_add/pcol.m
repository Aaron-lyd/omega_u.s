function [hf,hax,hp] = pcol(varargin)
% Plot a 2D matrix using pcolor. Returns handle to the figure.
%
%
% [hf,hax,hp] = pcol(ax, x, y, F, col_nan, ...)
% maps F using pcolor into axes ax, with x and y being the axes, and nan values
% coloured col_nan, and additional arguments passed directly to the construction of a new figure
% (only valid if ax is not provided).
% F can also be a string, in which case it is evaluated in the caller workspace, and title(F) is added.
% ax can also be 'current' to plot into the current axes (gca). 
% Outputs are handles to the figure, axes, and pcolor plot. 
% 
% Geoff Stanley. geoff.stanley@physics.ox.ac.uk
% SEE ALSO: figi

% Set defaults
x = [];
y = [];
NEWFIG = true;
col_nan = [1 1 1] * .2; % dark grey

narginchk(1,Inf);

% Check for x or 'current' argument
if isa(varargin{1},'matlab.graphics.axis.Axes')
  hax = varargin{1};
  hf = hax.Parent;
  varargin = varargin(2:end);
  NEWFIG = false;
elseif ischar(varargin{1}) && strncmpi(varargin{1}, 'current', 7)
  hax = gca();
  hf = gcf();
  varargin = varargin(2:end);
  NEWFIG = false;
end

% Check for x, y arguments
if length(varargin) >= 3 && isnumeric(varargin{1}) && isnumeric(varargin{2})
  % First two inputs are the grid for the x and y axes
  x = varargin{1};
  y = varargin{2};
  varargin = varargin(3:end);
end

% Extract field argument, evaluating it from caller workspace, if a string is passed.
fieldname = varargin{1};
varargin = varargin(2:end);

if ischar(fieldname)
  field = evalin('caller', fieldname);
else
  field = fieldname;
end
field = squeeze(field);
[ni,nj] = size(field);

% Transpose so first dim is horizontal (x), second dim is vertical (y)
field = [field' nan(nj,1); nan(1,ni+1)];

% Extract nan colour argument
if length(varargin) >= 1 && isnumeric(varargin{1}) && length(varargin{1}) == 3
  col_nan = varargin{1};
  varargin = varargin(2:end);
end

% Make new figure
if NEWFIG
  %aspect_ratio = 2.5;
  aspect_ratio = 5/3;
  h = 480;
  w = h * aspect_ratio;
  hf = figure('Position', [1920-w 1000-h w h], varargin{:});  % upper right corner of a 1920x1080 screen.
  hax = axes();
end
hold(hax, 'on'); % Do this before calling pcolor, or else it will destroy pre-existing colorbars.


% Pad coordinates by linear extrapolation, and pad data field with NaN's, to
% prep for pcolor with shading faceted / flat.
if isempty(x)
  x = 1 : (ni+1);
elseif isscalar(x)
  x = [x, x+1];
elseif isvector(x)
  x(end+1) = 2 * x(end) - x(end-1); % assume linear
else
  x = [x; 2*x(end,:) - x(end-1,:)];
  x = [x, 2*x(:,end) - x(:,end-1)];
end

if isempty(x)
y = 1 : (nj+1);
elseif isscalar(y)
  y = [y, y+1];
elseif isvector(y)
  y(end+1) = 2 * y(end) - y(end-1); % assume linear
else
  y = [y; 2*y(end,:) - y(end-1,:)];
  y = [y, 2*y(:,end) - y(:,end-1)];
end



hp = pcolor(hax, x, y, field);
hp.EdgeColor = 'none';

axis(hax, 'tight');
hax.Box = 'on';
hax.Color = col_nan; % nan colour.
%if isempty(findall(hf, 'Type', 'colorbar'))
%    % Add a colorbar if the figure doesn't already have one
%    colorbar(hax);
%end


hax.Layer = 'top';  % bring the axes (e.g. grid, and bounding box) to top.

if ischar(fieldname)
  title(fieldname, 'Interpreter', 'none');
end