function h = hist2D(varargin)
% hist2Dpatch is a modification of SMOOTHHIST2D to use patch objects to
% show the data, thereby allowing transparency and colours to be chosen,
% as well as to overlay multiple data sets with
% presumably different colours onto the same plot, by multiple calls to
% this command.
% Smoothing (lambda parameter) is removed.
% Some features are given through args, mostly allowing to pass plotting
% style information to 'patch'.
% Added a 3rd dimension to the plot: call this with additional parameters
%  ('3D', val), where `val` is the constant value in the extra Z dimension.
% You may also call this with ('3DY', val) to specify the value in the Y
%  dimension, and X(:,1) stays on the x axis, X(:,2) becomes the z axis.
% You may also call this with ('3DX', val) to specify the value in the X
%  dimension, and X(:,1) is moved to the y axis, X(:,2) to the z axis.
% SMOOTHHIST2D Plot a smoothed histogram of bivariate data.
%   SMOOTHHIST2D(X,LAMBDA,NBINS) plots a smoothed histogram of the bivariate
%   data in the N-by-2 matrix X.  Rows of X correspond to observations.  The
%   first column of X corresponds to the horizontal axis of the figure, the
%   second to the vertical. LAMBDA is a positive scalar smoothing parameter;
%   higher values lead to more smoothing, values close to zero lead to a plot
%   that is essentially just the raw data.  NBINS is a two-element vector
%   that determines the number of histogram bins in the horizontal and
%   vertical directions.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF) plots outliers in the data as points
%   overlaid on the smoothed histogram.  Outliers are defined as points in
%   regions where the smoothed density is less than (100*CUTOFF)% of the
%   maximum density.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,[],'surf') plots a smoothed histogram as a
%   surface plot.  SMOOTHHIST2D ignores the CUTOFF input in this case, and
%   the surface plot does not include outliers.
%
%   SMOOTHHIST2D(X,LAMBDA,NBINS,CUTOFF,'image') plots the histogram as an
%   image plot, the default.
%
%   Example:
%       X = [mvnrnd([0 5], [3 0; 0 3], 2000);
%            mvnrnd([0 8], [1 0; 0 5], 2000);
%            mvnrnd([3 5], [5 0; 0 1], 2000)];
%       smoothhist2D(X,5,[100, 100],.05);
%       smoothhist2D(X,5,[100, 100],[],'surf');
%
%   Reference:
%      Eilers, P.H.C. and Goeman, J.J (2004) "Enhancing scaterplots with
%      smoothed densities", Bioinformatics 20(5):623-628.

%   Copyright 2009 The MathWorks, Inc.
%   Revision: 1.0  Date: 2006/12/12
%
%   Requires MATLAB� R14.


[ax, args] = axescheck(varargin{:});
if isempty(ax)
    ax = axes();
end
hold(ax, 'on');
ax2 = [];

xdat = args{1};
ydat = args{2};
args = args(3:end);

nc = 256;

% If 'zerocolour' specified, use it, otherwise use default
ind = find(strncmpi(args, 'zerocol', 7));
if ind
    zerocolour = args(ind+1);    zerocolour = zerocolour{1};
else
    %zerocolour = [.8 .8 1]; %light blue
    zerocolour = [.8 .8 .8]; %grey
end


% If 'colormap' specified, use it, otherwise use @hot. For use only with 'plottype', 'image'.
ind = find(strcmpi(args, 'colormap'));
if ind
    cm = args(ind+1);    cm = cm{1};
    %assert(isa(cm, 'function_handle'), 'colormap value should be a function handle, e.g. @jet');
else
    cm = @hot;
end

if isa(cm, 'function_handle')
    cmap = feval(cm, nc);
else
    cmap = cm;
end

% If 'ROI' specified, use it, otherwise use full data
ind = find(strcmpi(args, 'ROI'));
if ind
    ROI = args(ind+1);    ROI = ROI{1};
    X = [flat(xdat(ROI)), flat(ydat(ROI))];
else
    X = [flat(xdat), flat(ydat)];
end
X(any(~isfinite(X),2),:) = [];


% If 'xlim' specified, use it, otherwise use default
ind = find(strcmpi(args, 'xlim'));
if ind
    xlim = args(ind+1);    xlim = xlim{1};
    X(X(:,1) < xlim(1), :) = [];  % Remove low  values of xdat from the data.
    X(X(:,1) > xlim(2), :) = [];  % Remove high values of xdat from the data.
end

% If 'ylim' specified, use it, otherwise use default
ind = find(strcmpi(args, 'ylim'));
if ind
    ylim = args(ind+1);    ylim = ylim{1};
    X(X(:,2) < ylim(1), :) = [];  % Remove low  values of ydat from the data.
    X(X(:,2) > ylim(2), :) = [];  % Remove high values of ydat from the data.
end

assert(size(X,1) > 0, 'You have managed to remove all your data!')

% If 'nbins' specified, use it, otherwise use default
ind = find(strcmpi(args, 'nbins'));
if ind
    nbins = args(ind+1);    nbins = nbins{1};
    if isscalar(nbins)
        nbins(2) = nbins(1);
    end
else
    nbins = [400 400];
end

% If 'OutlierShow' specified, use it, otherwise don't compute outliers any differently.
ind = find(strcmpi(args, 'Outlier'));
if ind
    outlier = args(ind+1);    outlier = outlier{1};
    default_outlier.frac = 0;
    default_outlier.size = 2;
    default_outlier.col = cmap(1,:);
    outlier = catstruct(default_outlier, outlier);
    %if ~isfield(outlier, 'ax')
    %    outlier.ax = axes();
    %end
else
    outlier.frac = 0;  % Show all data as normal, by default.
end

% If 'OutlierHide' specified, use it, otherwise don't compute outliers any differently.
ind = find(strcmpi(args, 'OutlierHide'));
if ind
    outlierhide = args(ind+1);    outlierhide = outlierhide{1};
else
    outlierhide = 0;  % Show all data, by default.
end

% If 'RescaleFcn' specified, use it. RescaleFcn is applied to relF, ie. values between 0 and 1 inclusive, and must be a vector function,
% ie. returning a matrix result given the matrix relF.
% Some useful options for rescaling function:
% @(x) x .^ 0.5 * 0.5    % scale up to better show the sparser data. Larger exponent = show sparse data more.
% @(x) (x > 0) .* 0.5    % show all data at half transparency.
% @(x) bsxfun(@rdivide, x.^1, max(x.^1, [], 1)); % Make each range of Y data have a bin with opacity 1.
ind = find(strcmpi(args, 'RescaleFcn'));
if ind
    rescalefcn = args(ind+1);    rescalefcn = rescalefcn{1};
else
    rescalefcn = @(x) x ;  % Identity scaling by default.
end

ind = find(strcmpi(args, 'EdgesX'));
if ind
    edges1 = args(ind+1);    edges1 = edges1{1};
    %ctrs1 = edges1(1:end-1) + .5*diff(edges1);
else
    minx = min(X(:,1));
    maxx = max(X(:,1));
    edges1 = linspace(minx, maxx, nbins(1)+1);
    %ctrs1 = edges1(1:end-1) + .5*diff(edges1);
end

ind = find(strcmpi(args, 'EdgesY'));
if ind
    edges2 = args(ind+1);    edges2 = edges2{1};
    %ctrs2 = edges2(1:end-1) + .5*diff(edges2);
else
    miny = min(X(:,2));
    maxy = max(X(:,2));
    edges2 = linspace(miny, maxy, nbins(2)+1);
    %ctrs2 = edges2(1:end-1) + .5*diff(edges2);
end

nbins = [length(edges1), length(edges2)] - 1;

[F,~,~,binx,biny] = histcounts2(X(:,1), X(:,2), [-Inf edges1(2:end-1) Inf], [-Inf edges2(2:end-1) Inf]);

n = size(X,1);
%bin = zeros(n,2);

% Reverse the columns of H to put the first column of X along the horizontal axis, the second along the vertical.
%[~,bin(:,2)] = histc(X(:,1), [-Inf edges1(2:end-1) Inf]);
%[~,bin(:,1)] = histc(X(:,2), [-Inf edges2(2:end-1) Inf]);
%F = accumarray(bin,1,nbins([2 1])); % was H = RHS.  I am doing no smoothing!
% F is nbins(1) by nbins(2) sized matrix.
% F(i,j) * n is the number of data points (from X) which are inside the (i,j) box.

% Eiler's 1D smooth, twice
%G = smooth1D(H,lambda);
%F = smooth1D(G',lambda)';
% % An alternative, using filter2.  However, lambda means totally different
% % things in this case: for smooth1D, it is a smoothness penalty parameter,
% % while for filter2D, it is a window halfwidth
%F = filter2D(F,lambda);

%relF = F./max(flat(F));
if outlier.frac > 0
    % Original approach: outliers = data from all those (2d) bins which have n*F =
    % number of data points in them being less than
    % outliershow * n*max(F(:)) = outliershow * (number of points in the most popular bin)
    %outliers = (relF(nbins(2)*(bin(:,2)-1)+bin(:,1)) < outliershow);
    
    % Geoff's alternative: outliers = data from all those bins with 1 data point,
    % then from all those bins with 2 data points, ... etc., up until the point
    % when we have collected (n*outliershow) outliers.
    % This way we can say easily, 'the plotted outliers are approximately the `outliershow` least
    % popular points.'
    % This is also better when n*max(F(:)) < 1/outliershow, in which case the above code
    % would pick out NO outliers.
    foo = F(sub2ind(size(F), binx, biny));
    for kk = 1 : max(flat(F))
        outliers = foo <= kk;
        if sum(outliers) > n * outlier.frac
            break
        end
    end
    
    % Hide bins where we will draw outliers:
    F(F <= kk) = 0;
end
%if outlierhide > 0
%    relF(relF < outlierhide) = 0;
%end

% Apply rescale function and ensure zero data gets the zerocolour while
% non-zero data gets colours from cm.
rF = rescalefcn(F).';
empty = (rF == 0);
rF(empty) = nan;

C = toTrueColor(rF, cmap);
rF(empty) = 0;
rF = cat(2, cat(1, rF, zeros(1,nbins(1))), zeros(nbins(2)+1,1));

h = surf(ax, edges1, edges2, rF, C);

% This style isn't as good. To make zero data invisible requires making rF == nan where it was zero
% but then this removes 4 faces around the vertex that was made to nan, not just one as it should!
%colormap(cmap);
%rF = cat(2, cat(1, rF, zeros(1,nbins(1))), zeros(nbins(2)+1,1));
%h = surf(ax, edges1, edges2, rF);

%h = surf(ax, edges1, edges2, rF, repmat(permute([0;0;0], [3 2 1]), [401 401 1]));
h.EdgeColor = 'none';
h.Parent.Color = zerocolour;
%hp.CDataMapping = 'direct';


% Set underlying axes to have same color limits and map as the histogram
% (in case user calls colorbar() afterwards)
ax.CLim = [min(rF(:)), max(rF(:))];
ax.Colormap = cmap;

% plot the outliers
if outlier.frac > 0
    %ax2 = axes('Position', ax.Position);
    %ax2.XLim = ax.XLim;
    %ax2.YLim = ax.YLim;
    %ax2.XTick = [];
    %ax2.YTick = [];
    %ax2.Color = 'none';
    hold(ax, 'on');
    plot(ax, X(outliers,1), X(outliers,2), '.', ...
        'MarkerEdgeColor', outlier.col, 'MarkerSize', outlier.size);
end


%     % plot a subsample of the data
%     Xsample = X(randsample(n,n/10),:);
%     plot(Xsample(:,1),Xsample(:,2),'bo');

%{
% For output: F = nbinsX by nbinsY
% F(i,j) gives number of data points in the bin centred at (ctrsX(i), ctrsY(j)).
F = F';

ind = find(strcmpi(args, 'overlay'));
if ind
    overlay = args(ind+1);      overlay = overlay{1};
    sum_over_2 = sum(F,2);
    signif = (sum_over_2 ./ sum(sum_over_2) > 0.001);
    if strcmpi(overlay, 'mean')
        % Weighted mean. F(i,j) where j is the dimension of xvar (e.g. density)
        % and i is the dimension of yvar (e.g. PV). So sum(F,2) sums F over the
        % y dimension, and (F * ctrsY') sums ctrsY weighted by the values of F
        % along the y dimension.
        F_bin_mean = (F * ctrs2') ./ sum_over_2 ;
        hp = plot(ctrs1(signif), F_bin_mean(signif), 'Color', facecol, 'LineWidth', 2);
    elseif strcmpi(overlay, 'mode')
        [~, I] = max(F,[],2);
        hp = plot(ctrs1(signif), ctrs2(I(signif)), 'Color', facecol, 'LineWidth', 2);
    end
    
    % Set legend to not display for the overlay line.
    if exist('hp', 'var')
        hAnnotation = get(hp,'Annotation');
        hLegendEntry = get(hAnnotation','LegendInformation');
        set(hLegendEntry,'IconDisplayStyle','off')
    end
    
    
end
if nargout > 1
    bin = bin(:,[2 1]);
end
%}
end

function C = toTrueColor(Z, cmap)
% Number of colors in colormap
N = size(cmap,1);

% Scale Z to range 1 .. N
minZ = min(Z(:));
maxZ = max(Z(:));
Z = (Z - minZ) * ((N-1) / (maxZ - minZ)) + 1;

% Interpolate colormap:
C = interp1((1:N)', cmap, Z, 'linear');
end

%--------------------------------------------------------------------------
%{
function Z = smooth1D(Y,lambda)
[m,n] = size(Y);
E = eye(m);
D1 = diff(E,1);
D2 = diff(D1,1);
P = lambda.^2 .* D2'*D2 + 2.*lambda .* D1'*D1;
Z = (E + P) \ Y;
% This is a better solution, but takes a bit longer for n and m large
% opts.RECT = true;
% D1 = [diff(E,1); zeros(1,n)];
% D2 = [diff(D1,1); zeros(1,n)];
% Z = linsolve([E; 2.*sqrt(lambda).*D1; lambda.*D2],[Y; zeros(2*m,n)],opts);
%}

%--------------------------------------------------------------------------
%{
function Z = filter2D(Y,bw)
z = -1:(1/bw):1;
k = .75 * (1 - z.^2); % epanechnikov-like weights
k = k ./ sum(k);
Z = filter2(k'*k,Y);
%}

function out = flat(in)
out = in(:);
end