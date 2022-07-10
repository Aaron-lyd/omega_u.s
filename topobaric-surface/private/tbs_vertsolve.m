function [p, s, t] = tbs_vertsolve(Sppc, Tppc, P, BotK, s, t, p, branchmap, d_fn, s0, t0, tolp, DP) %#codegen
%TBS_VERTSOLVE  Root finding of pressure or depth that matches equation of
%               state with multivalued function.
%
%
% [p, s, t] = tbs_vertsolve(Sppc, Tppc, P, BotK, s, t, p, branch, d_fn, s0, t0, tolp, DP)
% finds the pressure or depth p (within tolerance tolp) and its associated
% practical / Absolute salinity s and potential / Conservative temperature
% t, of a surface on which delta equals that determined by the multivalued
% function d_fn, in an ocean whose practical / Absolute salinity and
% potential / Conservative temperature as functions of pressure or depth P
% are given by piecewise polynomials whose coefficients are Sppc and Tppc,
% and whose knots are P.  The number of valid data points in each water
% column is given by BotK.  The equation of state is given by eos.m in the
% path, taking S, T, P as its 3 inputs. delta is the in-situ density
% anomaly or specific volume anomaly, defined as eos(s,t,p) - eos(s0,t0,p)
% where s,t are S,T interpolated from P to p, and s0, t0 are reference
% values.  In each water column I, the branch of the multivalued function
% for delta as a function of p is given by d_fn(:,b) where b =
% branchmap(I).  This branch is a polynomial of p between d_fn(1,b) and
% d_fn(2,b) with coefficients given by d_fn(3:end,b), and a linear
% extension of this polynomial outside this domain. Solutions are sought in
% the interval [A,B] where A = max(P(1,I), d_fn(1,b) - DP) and B =
% min(P(BotK(I),I), d_fn(2,b) + DP), i.e. within the valid P range of the
% local water column AND within DP of the polynomial domain of the local
% branch. (If P is a vector, the I indexing of P is dropped.)  Limiting
% this search domain with DP helps prevent the surface from jumping between
% multiple solutions in weakly stratified waters such as around Antarctica.
% The inputs s and t are not used, but provided so these variables may be
% manipulated in-place.
%
%
% --- Input:
% Sppc [O, K-1, N]: coefficients for piecewise polynomial for practical
%                   / Absolute Salinity in terms of P
% Tppc [O, K-1, N]: coefficients for piecewise polynomial for potential
%                   / Conservative Temperature in terms of P
% P [K, N]: knots for the pressure [dbar] or depth [m] of the casts
% BotK [1, N]: number of valid data points on each cast
% s [1, N]: initial practical / Absolute salinity on the initial surface
% t [1, N]: initial potential / Conservative temperature on the initial surface
% p [1, N]: initial pressure [dbar] or depth [m] on the initial surface
% branchmap [1, N]: branch that each cast belongs to.
%                   Must have 1 <= branch(i) <= A, for all 1 <= i <= N.
% d_fn [D+3, A]: domain and coefficients for each branch's polynomial
%             for delta as a function of p, the b'th branch being defined
%             by d_fn(:,b), which can be evaluated by pvaln or pvallin.
% s0 [1, 1]: reference S value for delta
% t0 [1, 1]: reference T value for delta
% tolp [1, 1]: tolerance on pressure [dbar] or depth [m] for root finding solver
% DP [1, 1]: seek solutions in the domain of the local d_fn branch expanded
%            by DP in both directions.
%
% Note: O is the order of the piecewise polynomials down each cast
%       K is the maximum number of knots in these piecewise polynomials,
%           i.e. the maximum number of bottles in any cast
%       N is the number of water columns (possibly including land).
%       A is the number of branches for the multivalued function d_fn.
%       D is the degree of each branch's polynomial
%
% Note: P must increase along its first dimension.
%
% Note: P can have size [K, 1], in which case it is used for each cast.
%
% Note: variables can actually be higher dimensional, e.g. N = [ni, nj],
%       and p can be any dimensional matrix, so long as it has N elements
%       in total.
%
% Note: BotK should be given by
%           BotK = squeeze(sum(isfinite(S), 1));
%
%
% --- Output:
% p [same as input p]: pressure or depth of the updated surface
% s [same as input p]: practical / Absolute salinity of the updated surface
% t [same as input p]: potential / Conservative temperature of the updated surface

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


N = numel(p);
Pmat = ~isvector(P);

% Loop over each cast
for n = 1:N
  b = branchmap(n);
  k = BotK(n);
  if b > 0 && k > 1
    
    % Select this water column
    Sppcn = Sppc(:,1:k-1,n);
    Tppcn = Tppc(:,1:k-1,n);
    if Pmat
      Pn = P(1:k,n);
    else
      Pn = P((1:k).'); % .' is for codegen, so P and (1:k).' both column vectors
    end
    
    lb = max(Pn(1), d_fn(1,b) - DP);
    ub = min(Pn(k), d_fn(2,b) + DP);
    
    % Search for a sign-change, expanding outward from an initial guess
    [lb, ub] = fzero_guess_to_bounds(@myfcn, p(n), lb, ub, ...
      Sppcn, Tppcn, Pn, d_fn(:,b), s0, t0);
    
    if ~isnan(lb)
      % A sign change was discovered, so a root exists in the interval.
      % Solve the nonlinear root-finding problem using Brent's method
      p(n) = fzero_brent(@myfcn, lb, ub, tolp, ...
        Sppcn, Tppcn, Pn, d_fn(:,b), s0, t0);
      
      % Interpolate S and T onto the updated surface
      [s(n),t(n)] = ppc_val2(Pn, Sppcn, Tppcn, p(n));
    else
      p(n) = nan;
      s(n) = nan;
      t(n) = nan;
    end
    
  else
    % This will ensure s,t,p all have the same nan structure
    p(n) = nan;
    s(n) = nan;
    t(n) = nan;
  end
  
end

end


function out = myfcn(p, Sppc, Tppc, P, d_fn, s0, t0)
% The difference between delta evaluated (a) using the local branch of the
% multivalued function, and (b) using the equation of state with the local
% water properties.

% Evaluate water properties on the surface
[s,t] = ppc_val2(P, Sppc, Tppc, p);

% Evaluate the delta difference
out = pvallin(d_fn, p) - ( eos(s, t, p) - eos(s0, t0, p) );

end
