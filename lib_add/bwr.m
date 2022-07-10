function varargout = bwr(c,ncol)
% Set the colour axis to blue white red with magnitude c
% Output is just to assist using this in anonymous functions, eg
% plt = @(fn,a) {subplot(3,3,a); figi(gU.(fn)(:,:,k), 'currentfig'); bwr(nanrms(flat(gU.(fn)(:,:,k)))); title(fn)};
if nargin < 2 || isempty(ncol)
    ncol = 128;
end
caxis([-1 1]*c);
colormap(gca, bluewhitered(ncol));
varargout = cell(nargout,1);