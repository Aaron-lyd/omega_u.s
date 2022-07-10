function [hax, h] = pcolormask(varargin)
% PCOLOR_NAN_MASK: make a pcolor plot with specified colours for nan mask.
% INPUT:
% X [vector length M]: the horizontal axis
% Y [vector length N]: the vertical axis
% mask [matrix M by N]: Where mask ~= 0, the pixel will be coloured with maskcol.
% maskcol: a RGB triplet (3 element vector)
% Can also pass an axis handle as a the first argument.
%
% OUTPUT:
% ax = handle to the axes
% h = handle to the surf
%
% Geoff Stanley
% geoff.stanley@physics.ox.ac.uk



[hax,args,nargs] = axescheck(varargin{:});
if nargs < 4
    error(message('MATLAB:narginchk:notEnoughInputs'));
end

X = args{1};
Y = args{2};
mask = args{3};
maskcol = args{4};
if nargs < 5
    z = 0;
else
    z = args{5};
end

if isempty(hax)
    figure;
end
hax = newplot(hax);

[nx, ny] = size(mask);
if isempty(X)
    X = 1:nx;
end
if isempty(Y)
    Y = 1:ny;
end
X = X(:);
Y = Y(:);

mask = mask.';
%z = [double(mask ~= 0), zeros(nx,1); zeros(1,ny+1)];
%z(z == 0) = nan;
C = nan(nx*ny,3);
C(mask,1) = maskcol(1);
C(mask,2) = maskcol(2);
C(mask,3) = maskcol(3);
C = reshape(C, [ny nx 3]);
C = [C, zeros(ny,1,3); zeros(1,nx+1,3)];
h = surf(hax, [X; Inf], [Y; Inf], z + zeros(ny+1,nx+1), C);
h.EdgeColor = 'none';
view(hax,0,90);
