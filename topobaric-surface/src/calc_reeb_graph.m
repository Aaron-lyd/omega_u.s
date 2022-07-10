function [p, varargout] = calc_reeb_graph(p, WRAP, OPTS)
%CALC_REEB_GRAPH  Calculate the Reeb graph of p
%
%
% [p, RG] = calc_reeb_graph(p, WRAP, OPTS)
% performs the following steps:
% 1, fill small holes with high data in the given rectilinear data p (if
%    OPTS.FILL_PIX is positive or OPTS.FILL_IJ is not empty);
% 2, build a simplical mesh given rectilinear data in p;
% 3, perturb the data p on that mesh until all values are unique;
% 4, calculate the Reeb graph of p on that mesh, using the ReCon software
%    (Doraiswamy and Natarajan, 2013);
% 5, post-process the Reeb graph to take out any holes that were filled with
%    high data
% 6, simplify the Reeb graph (if it has more arcs than
%    OPTS.SIMPLIFY_ARC_REMAIN).
%
%
% --- Input:
% p [ni, nj]: rectilinear data, containing exactly ONE connected component
%             (where each pixel has 4 neighbours, not 8), as can be found
%             using bfs_conncomp(). 
% WRAP [2 element array]: determines which dimensions are treated periodic
%                         [logical].  Set WRAP(i) to true when periodic in 
%                         the i'th lateral dimension(i=1,2).
% OPTS [struct]: options. See topobaric_surface documentation for details.
%
%
% --- Output:
% p [ni, nj]: identical to the input p but possibly having some values
%             perturbed by an amount near machine precision.  The input
%             and output p have identical finite vs NaN structure. 
% RG [struct]: The Reeb graph. See topobaric_surface.m for details.
% timer_recon:  time from tic/toc spent running ReCon. 
%
% If more than three outputs are requested, the first output is p and the
% rest of the outputs are the following fields of RG: arc_from, arc_to,
% arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs,
% nNodes, wet, n_casts.  Finally, timer_recon is the last output. 
%
%
% --- References:
% Doraiswamy, H. & Natarajan, V. Computing Reeb Graphs as a Union of
% Contour Trees. IEEE Transactions on Visualization and Computer Graphics
% 19, 249?262 (2013).

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com



[ni,nj] = size(p);

% wet: a logical map where p is valid. As this function proceeds, wet ==
% isfinite(p) is maintained.  wet at the end of this function is identical
% to wet at the start of this function.  That is, isfinite(output p) ==
% isfinite(input p).
wet = isfinite(p);  


OPTS = catstruct(tbs_defaults(), OPTS);


% --- Pre-processing 1: Fill in small or specific holes
DO_FILL = OPTS.FILL_PIX > 0 || ~isempty(OPTS.FILL_IJ);
if DO_FILL
    
    Ground = ~wet;
    
    % Calculate connected components of Ground.  Use 8 connectivity here!
    neigh = grid_adjacency([ni,nj], 8, WRAP);
    [qu,qts,~,~,L] = bfs_conncomp(Ground, neigh, []);
    
    % Set to false (will not fill) for land regions that have > specified number of pixels:
    if OPTS.FILL_PIX > 0
        npix = diff(qts);
        for c = find(npix(:).' > OPTS.FILL_PIX)
            Ground(qu(qts(c) : qts(c+1)-1)) = false;
        end
    end
    
    % Set to false (will not fill) for any land regions that are connected
    % to specified locations (e.g. the continents and major islands)
    if ~isempty(OPTS.FILL_IJ)
        
        for I = 1:size(OPTS.FILL_IJ,1)
            i = OPTS.FILL_IJ(I,1);
            j = OPTS.FILL_IJ(I,2);
            if Ground(i,j)
                c = L(i,j);
                Ground(qu(qts(c) : qts(c+1)-1)) = false ;
            end
        end
    end
    
    % Now fill in: It won't matter what values are put in there, so long as
    % they are all greater than all surrounding pixels.
    fill = find(Ground);
    offset = max(p(:));
    p(fill) = (1+offset) : (length(fill) + offset);
    
    % And update wet, to maintain wet == isfinite(p) 
    wet(fill) = true;
end


% --- Pre-processing 2: Make simplical decomposition
[verts,faces,~,nV] = latlon2vertsfaces(p,WRAP,[],[],OPTS.DECOMP,false);

% Last argument (select) could be OPTS.REF_IJ as below, but we use false,
% since we've pre-selected one connected region.
%[verts,faces,~,nV] = latlon2vertsfaces(p,WRAP,[],[],OPTS.DECOMP,OPTS.REF_IJ);

% Ensure all vertices have unique values
[verts(:,3), perturbed] = perturb_until_unique(verts(:,3));

mytic = tic(); % DEV
% --- Processing: the Reeb Graph (RG)
% Note: Recon (version 1.0) can only handle surfaces with one connected
% component.
RAA = vgl.iisc.recon.incore.ReconAlgorithmAug;
faces = int32(faces-1); % -1 because Java wants 0-indexing, but MATLAB does 1-indexing
RAA.run(verts(:,1), verts(:,2), zeros(size(verts,1),1), verts(:,3), ...
    faces(:,1), faces(:,2), faces(:,3));
clear faces

% Bring the RG into MATLAB as a struct.
% First argument (nV): causes vertices numbered above nV to be left out of arc_segment
% These are not only the extra vertices formed from the simplical decomposition, but also
% extra vertices added inside ReCon when splitting the input into loop free regions.
% Second argument (1): shift everything by +1, for MATLAB 1 indexing.
RG = struct(RAA.pack(nV,1));
clear RAA

timer_recon = toc(mytic);

% Handle the special case of only one arc.
if RG.nArcs == 1
  RG.arc_segment = {RG.arc_segment};
end

% Dereference RG for faster execution. Also convert to double (which is not slower overall)
arc_from    = double(RG.arc_from);
arc_to      = double(RG.arc_to);
arc_segment = cellfun(@double, RG.arc_segment, 'UniformOutput', false);
node_fn     = RG.node_fn;
node_next   = cellfun(@double, RG.node_next, 'UniformOutput', false);
node_prev   = cellfun(@double, RG.node_prev, 'UniformOutput', false);
node_type   = RG.node_type; % Leave as int8
node_v      = double(RG.node_v);
nArcs       = RG.nArcs;
nNodes      = RG.nNodes;


% --- Post-processing 1: Update p to match the verts used in the Reeb graph
% (if they were perturbed)
if perturbed
    % Just update those vertices that were from the original rectangular grid
    p(wet) = verts(1:nV,3);
end
clear verts


% --- Post-processing 2: if holes were filled, fix node_v and arc_segments
if DO_FILL
    % Adjust node_v
    keep_verts = ~Ground(wet);
    vertmap = cumsum(keep_verts);
    vertmap(~keep_verts) = 0; % Mark, as 0, any vertices (and nodes via node_v) that were inside filled regions
    if lower(OPTS.DECOMP(1)) == 'd'
        node_v = vertmap(node_v);
    else
        node_v(node_v <= nV) = vertmap(node_v(node_v <= nV));
    end
    
    % Remove missing vertices from arc_segment, and update the vertex indices
    for e = 1:nArcs
        arc_segment{e} = vertmap(arc_segment{e}(keep_verts(arc_segment{e})));
    end
    
    
    % Finally, reset any pixels that were filled to nan
    p(fill) = nan;
    
    % And update wet, to maintain wet == isfinite(p) 
    wet(fill) = false;
    
    % And update the number of casts
    nV = nV - length(fill);
    
end

% --- Post-processing 3: Simplify the Reeb Graph:
if OPTS.SIMPLIFY_ARC_REMAIN < nArcs
    assert(OPTS.DECOMP(1) == 'd', 'Reeb graph simplification with CROSS simplical decomposition disabled')
    [arc_from, arc_to, arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs, nNodes] ...
        = simplify_reeb_graph( ...
        arc_from, arc_to, arc_segment, node_prev, node_next, node_type, node_fn, node_v, nArcs, ...
        OPTS.SIMPLIFY_ARC_REMAIN, OPTS.SIMPLIFY_WEIGHT_PERSIST ...
        );
end

% --- Integrity checks:
%{
assert(all(wet(:) == isfinite(p(:))), 'wet differs from where p is finite');

livenodes = 0 < node_v & node_v <= nV;

p_ = p(wet);
chk = p_(node_v(livenodes)) == node_fn(livenodes);
assert(all(chk), 'Check failed: node_fn, node_v, and p_ are mismatched.');
clear livenodes chksimplify_reeb_graph

for e = 1:nArcs
    assert( sum(e == node_next{arc_from(e)}) == 1, 'Arc from mismatch!' )
    assert( sum(e == node_prev{arc_to(  e)}) == 1, 'Arc to   mismatch!' )
end
for n = 1:nNodes
    for e = node_next{n}.'
        assert( arc_from(e) == n, 'node next broken!');
    end
    for e = node_prev{n}.'
        assert( arc_to(e) == n, 'node prev broken!');
    end
end
%}

if nargout <= 3
    % Pack output into a struct
    RG = struct();
    RG.arc_from    = arc_from;
    RG.arc_to      = arc_to;
    RG.arc_segment = arc_segment;
    RG.node_fn     = node_fn;
    RG.node_next   = node_next;
    RG.node_prev   = node_prev;
    RG.node_type   = node_type;
    RG.node_v      = node_v;
    RG.nArcs       = nArcs;
    RG.nNodes      = nNodes;
    RG.wet         = wet;
    RG.n_casts     = nV;
    varargout{1} = RG;
    varargout{2} = timer_recon;
else
    % Return a long list of outputs
    varargout{1} = arc_from;
    varargout{2} = arc_to;
    varargout{3} = arc_segment;
    varargout{4} = node_prev;
    varargout{5} = node_next;
    varargout{6} = node_type;
    varargout{7} = node_fn;
    varargout{8} = node_v;
    varargout{9} = nArcs;
    varargout{10} = nNodes;
    varargout{11} = wet;
    varargout{12} = nV;
    varargout{13} = timer_recon;
end

