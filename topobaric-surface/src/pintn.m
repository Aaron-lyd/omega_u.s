function F = pintn(f,a,b)
%PINTN  Integrate many polynomials.
%
%
% F = pintn(f)
% evaluates the indefinite integral of the polynomials f, with lower
% integration limit given by the lower limit of the domain of f, producing
% new polynomials F (one order higher than f and valid on the same domain
% as f). Evaluating F at x (using pvaln) gives the integral of f from the
% lower limit of the domain to x. Each column of f and F represents a
% polynomial.
%
% F = pintn(f,a)
% as above, but the lower limit of the integral is a, which is a scalar or
% a vector with as many elements as there are columns in f.
%
% F = pintn(f,a,b)
% evaluates the definite integral of f from a to b. This simply
% evaluates F = pintn(f,a) at b.
%
% Though the domain of f is specified, this is ignored for execution speed.
%
% Any of f, a, or b can have a singleton second dimension, in which case
% that polynomial or integration limit is used for all integrals.
%
% There are three input/output cases.
%
%
% --- Input (a):
% f [N+2, L]: L polynomials of order N (may have L = 1)
%
%
% --- Output (a):
% F [N+2, L]: F(:,l) is the indefinite integral of f(:,l) from f(1,l)
%
%
% --- Input (b):
% f [N+2, L]: L polynomials of order N (may have L = 1)
% a [1  , L]: lower limits of integration (may have L = 1)
%
%
% --- Output (b):
% F [N+2, L]: F(:,l) is the indefinite integral of f(:,l) from a(1,l)
%
%
% --- Input (c):
% f [N+2, L]: L polynomials of order N (may have L = 1)
% a [M, L]: lower limits of integration (may have L = 1)
% b [M, L]: lower limits of integration (may have L = 1)
%
%
% --- Output (c):
% F [M, L]: F(m,l) is the definite integral of f(:,l) from a(m,l) to b(m,l)

% Author(s) : Geoff Stanley
% Email     : g.stanley@unsw.edu.au
% Email     : geoffstanley@gmail.com


% Np2-2 is the polynomial's order. (e.g. Np2==4 for linear)
[Np2, M] = size(f);

% F is one order higher than f, so add a row to F
F = vertcat(f, zeros(1,M));

% Integration from lower limit of domain
F(3:Np2-1,:) = F(3:Np2-1,:) ./ (Np2-2:-1:2).' ;

if nargin >= 2
    % Make indefinite integral from lower value of a
    % (forgo check that a is in the domain)
    F(Np2+1,:) = -pvaln(F, a);
    
    if nargin == 3
        % Definite integral from a to b
        % (forgo check that b is in the domain)
        F = pvaln(F, b);
    end
    
end