function output = div_dot(ax, ay, bx, by, A_x, A_y, A)

im1 = @(F) circshift(F, [+1 0]);
jm1 = @(F) circshift(F, [0 +1]);
ip1 = @(F) circshift(F, [-1 0]);
jp1 = @(F) circshift(F, [0 -1]);

cx = ax .* bx .* (A_x);
cy = ay .* by .* (A_y);
output = (cx + ip1(cx) + cy + jp1(cy)) ./ (2 * A);