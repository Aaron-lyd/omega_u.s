function [y,yx] = pchipdqn(x,X,Y) %#codegen
%PCHIPDQN  Piecewise Cubic Hermite Interpolation and its Derivative along all columns.
%
%
% [y,yx] = pchipdqn(x,X,Y)
% uses a PCHIP to interpolate each column of Y as a function of each column
% of X, at each column of x. That is, y(l,m) interpolates Y(:,m), as a
% function of X(:,m), at U(l,m), and yx(l,m) is the derivative of this
% function. No extrapolation is done: when x(l,n) is out of bounds of
% X(:,n), y(l,n) and yx(l,n) are NaN.
%
%
% --- Input:
% x [L x N], the interpolant evaluation sites
% X [M x N], the data sites for the functions
% Y [M x N], the values of the first function at the data sites
%
%
% --- Output:
% y [L x N] , the interpolant values for the first function
% yx [L x N], the derivative values for the first function
%
%
% --- Notes:
% X(:,n) must be monotonically increasing for all n.
%
% NaN's in X are treated as +Inf (and as such must come at the end of each
% column).
%
% Any input can be higher-dimensional, so long as the product of its
% dimensions excluding the first dimension is N.
%
% Any input can have a singleton second dimension (N = 1), in which case
% that single value is used for each interpolation problem.
%
% Even if L == 1, x does need a leading singleton dimension.
%
% If L == 1, this dimension is squeeze()'d out of y.
%
%
% --- Code generation:
% PCHIPDQN can be compiled into a MEX executible using codegen, as follows:
% >> codegen('pchipdqn', '-args', {zeros(L,N), zeros(M,N), zeros(M,N)}, '-o', 'pchipdqn_mex');
% To call the mex function rather than this, call pchipdqn_mex().
% A more flexible example, using an upper bound for certain dimensions of
% input, is as follows:
% >> type_x = coder.typeof(0, [L P Q], [false true true]);
% >> type_X = coder.typeof(0, [M P Q], [false true true]);
% >> codegen('pchipdqn', '-args', {type_x, type_X, type_X}, '-o', 'pchipdqn_mex');
%
%
% Author    : Geoff Stanley
% Email     : geoffstanley@gmail.com
% Version   : 1.1
% History   : 05/09/2019 - initial release
%           : 25/10/2019 - minor speed and documentation updates

szX = size(X);
szY = size(Y);
XN = prod(szX(2:end));
YN = prod(szY(2:end));

% If x given as [N], add leading singleton dimension.
% This has been removed, as a necessity for codegen.
szx = size(x);
% xn = numel(x);
% if szx(1) ~= 1 && xn == max(XN, YN)
%     x = reshape(x, [1, szx]);
%     szx = size(x);
% end
xN = prod(szx(2:end));
M = szX(1);
L = szx(1);

[N,ii] = max([xN, XN, YN]);
assert(szY(1) == M, 'X and Y must have the same number of rows.');
assert(xN == 1 || xN == N, 'x must have N=%d columns or 1 column', N);
assert(XN == 1 || XN == N, 'X must have N=%d columns or 1 column', N);
assert(YN == 1 || YN == N, 'Y must have N=%d columns or 1 column', N);

szy = {szx(2:end), szX(2:end), szY(2:end)};
szy = [L, szy{ii}];
y = nan(szy);
yx = nan(szy);

% Variables to index the 1'st or n'th column:
x01 = double(xN > 1);
MX = M * double(XN > 1);
MY = M * double(YN > 1);
nX = 0 ;
nY = 0 ;

% Pre-assign sizes for PCHIP variables.
h = zeros(3,1);
deltaY = zeros(3,1);
dY = zeros(2,1);

% Loop over columns, each of which is an interpolation problem
for n = 0:N-1
    if isnan(Y(nY + 2)) || isnan(X(nX + 2))
        nX = nX + MX;
        nY = nY + MY;
        continue
    end
    for l = 1:L
        xln = x(n*L*x01 + l); % x(l,n) or x(l,1) as appropriate
        
        if isnan(xln) || xln < X(nX + 1) || xln > X(nX + M)
            
            y(n*L + l) = nan;
            
        else
            
            % Binary search, to find i for which X(i-1,n) < xln <= X(i,n)
            i = 2; % Result will be >= 2 always
            k = M; % Result will be <= M always
            while (i < k)
                j = floor((i + k) / 2);
                if X(nX + j) < xln % [X(j,n) or X(j,1) as appropriate]  < xln
                    i = j + 1;
                else
                    k = j;
                end
            end
            
            % Check whether xln is adjacent to the start or end of this X
            at_start = (i == 2);
            at_end = (i == M) || isnan(X(nX + i+1)) || isnan(Y(nY + i+1));
            
            if at_start && at_end
                % . X(1) < x < X(2) .   Revert to Linear Interpolation
                d = (xln - X(nX + i-1)) ./ (X(nX + i) - X(nX + i-1));
                y(n*L + l) = Y(nY+i-1) * (1-d) + Y(nY+i) * d;
                
            else
                if at_start
                    % . X(1) < x <= X(2) < X(3) ...
                    h(2) = X(nX+i  ) - X(nX+i-1);
                    h(3) = X(nX+i+1) - X(nX+i  );
                    deltaY(2) = (Y(nY+i  ) - Y(nY+i-1)) / h(2);
                    deltaY(3) = (Y(nY+i+1) - Y(nY+i  )) / h(3);
                    
                    %  Noncentered, shape-preserving, three-point formula:
                    dY(1) = ((2*h(2) + h(3)) * deltaY(2) - h(2) * deltaY(3)) / (h(2)+h(3));
                    if sign(dY(1)) ~= sign(deltaY(2))
                        dY(1) = 0;
                    elseif (sign(deltaY(2)) ~= sign(deltaY(3))) && (abs(dY(1)) > abs(3*deltaY(2)))
                        dY(1) = 3 * deltaY(2);
                    end
                    
                    % Standard PCHIP formula
                    if sign(deltaY(2)) * sign(deltaY(3)) > 0
                        w1 = 2*h(3) +   h(2);
                        w2 =   h(3) + 2*h(2);
                        dY(2) = (w1+w2) / (w1/deltaY(2) + w2/deltaY(3));
                    else
                        dY(2) = 0;
                    end
                    
                elseif at_end
                    % ... X(i-2) < X(i-1) < x <= X(i) .
                    h(1) = X(nX+i-1) - X(nX+i-2);
                    h(2) = X(nX+i  ) - X(nX+i-1);
                    deltaY(1) = (Y(nY+i-1) - Y(nY+i-2)) / h(1);
                    deltaY(2) = (Y(nY+i  ) - Y(nY+i-1)) / h(2);
                    
                    % Standard PCHIP formula
                    if sign(deltaY(1)) * sign(deltaY(2)) > 0
                        w1 = 2*h(2) +   h(1);
                        w2 =   h(2) + 2*h(1);
                        dY(1) = (w1+w2) / (w1/deltaY(1) + w2/deltaY(2));
                    else
                        dY(1) = 0;
                    end
                    
                    %  Noncentered, shape-preserving, three-point formula:
                    dY(2) = ((h(1) + 2*h(2)) * deltaY(2) - h(2) * deltaY(1)) / (h(1)+h(2));
                    if sign(dY(2)) ~= sign(deltaY(2))
                        dY(2) = 0;
                    elseif (sign(deltaY(2)) ~= sign(deltaY(1))) && (abs(dY(2)) > abs(3*deltaY(2)))
                        dY(2) = 3 * deltaY(2);
                    end
                    
                else
                    % ... X(i-2) < X(i-1) < x <= X(i) < X(i+1) ...
                    h(1) = X(nX+i-1) - X(nX+i-2); % Way faster to do this
                    h(2) = X(nX+i  ) - X(nX+i-1); % than  
                    h(3) = X(nX+i+1) - X(nX+i  ); % diff(X(nx+i-2:nX+i+1))
                    deltaY(1) = (Y(nY+i-1) - Y(nY+i-2)) / h(1);
                    deltaY(2) = (Y(nY+i  ) - Y(nY+i-1)) / h(2);
                    deltaY(3) = (Y(nY+i+1) - Y(nY+i  )) / h(3);
                    
                    % Standard PCHIP formula
                    for j = 1:2
                        if sign(deltaY(j)) * sign(deltaY(j+1)) > 0
                            w1 = 2*h(j+1) +   h(j);
                            w2 =   h(j+1) + 2*h(j);
                            dY(j) = (w1+w2) / (w1/deltaY(j) + w2/deltaY(j+1));
                        else
                            dY(j) = 0;
                        end
                    end
                end
                
                % Polynomial coefficients for this piece
                cY = (3*deltaY(2) - 2*dY(1) - dY(2)) / h(2);
                bY = (dY(1) - 2*deltaY(2) + dY(2)) / h(2)^2;
                
                % Evaluate cubic interpolant
                s = xln - X(nX + i-1);
                y(n*L + l) = Y(nY + i-1) + s * (dY(1) + s * (cY + s * bY));
                
                % Derivatives can be evaluated easily:
                yx(n*L + l) = dY(1) + s*(2*cY + 3*s*bY); % first derivative
                % y2(n*L + l) = 2*cY + 6*s*bY;            % second derivative
                % y3(n*L + l) = 6*bY;                      % third derivative
            end
        end
    end % for l
    nX = nX + MX;
    nY = nY + MY;
end

% Remove leading dimensions if L == 1, but leave row vectors as row vectors
y = squeeze(y);
yx = squeeze(yx);
