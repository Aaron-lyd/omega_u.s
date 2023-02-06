function sol = nonzero_mean_brent_lsqlin(mat, rhs, mr, pinval)
% Find the least-squares solution `sol` to the possibly overdetermined
% matrix problem `mat * sol = rhs` subject to the constraint that
% `sum(sol)` equals a certain value which happens to yield `sol(mr) =
% pinval`. This entails a root finding problem on a nonlinear function F(c)
% = sol_mr, which is the `mr`'th component of sol, and sol is the
% least-sqares solution to `mat * sol = rhs` subject to `sum(sol) = c`.

% Figure out what mean value will give Phi = 0 at ref cast, by Brent's method withLSQlin
opts = optimoptions('lsqlin');
opts.Display = 'off';
% opts.StepTolerance = 1e-15;
% opts.OptimalityTolerance = 1e-15;
% opts.ConstraintTolerance = 1e-15;

% Aaron OPTS
opts.StepTolerance = 1e-12;
opts.OptimalityTolerance = 2e-14;
opts.ConstraintTolerance = 1e-8;

N = size(mat, 2);
Aeq = ones(1,N);  % (later, use area weighting)
fac = sum(Aeq);

a = pinval * fac; % Initial guess for a root
f1 = F(a, mat, rhs, Aeq, mr, pinval, opts);

% Second guess for a root
b = -f1 * fac * 1.2;  % Try to overshoot the zero by a little bit (~20%)
f2 = F(b, mat, rhs, Aeq, mr, pinval, opts);
if b < a  % swap
  tmp = a;
  a = b;
  b = tmp;
end

% Check for a sign change between a and b.  If not, move a left and b right, and test again.
if f1 * f2 >= 0  % no sign change
  for rep = 1:100
    gap = b - a;
    a = a - gap * 0.2;
    b = b + gap * 0.2;
    if F(a, mat, rhs, Aeq, mr, pinval, opts) * F(b, mat, rhs, Aeq, mr, pinval, opts) < 0
      break
    end
  end
end

% Now solve the nonlinear root finding problem, to get c at which myfun_lsqlin(c, ...) = 0.
tol = 1e-16;
c = fzero_brent(@F, a, b, tol, mat, rhs, Aeq, mr, pinval, opts);

% Get the full solution of the matrix problem, with mean = c.
%[sol, ~, ~, flag] = lsqlin(mat, rhs, [], [], Aeq, c, [], [], [], opts);
[~, sol] = F(c, mat, rhs, Aeq, mr, pinval, opts);

%assert(sol(mr) <= tol);

end

function [f, sol] = F(c, mat, rhs, Aeq, mr, pinval, opts)
[sol, ~, ~, flag] = lsqlin(mat, rhs, [], [], Aeq, c, [], [], [], opts);
if flag ~= 1; warning('flag!'); end
f = sol(mr) - pinval;
end