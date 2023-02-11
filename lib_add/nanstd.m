function y = nanstd(x)

x = x(isfinite(x));
x = x - mean(x);
y = sqrt(mean(x .* conj(x)));
