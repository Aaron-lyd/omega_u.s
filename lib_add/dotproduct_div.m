function [c_A, c_A2] = dotproduct_div(c, dx, dy, A_x, A_y, A)
im1 = @(F) circshift(F, [+1 0]);
jm1 = @(F) circshift(F, [0 +1]);
ip1 = @(F) circshift(F, [-1 0]);
jp1 = @(F) circshift(F, [0 -1]);

cx = (c - im1(c))./dx;
cy = (c - jm1(c))./dy;

cx2 = cx .* cx .* (A_x);
cy2 = cy .* cy .* (A_y);
c_A2 = (cx2 + ip1(cx2) + cy2 + jp1(cy2)) ./ (2 * A);
c_A = sqrt(c_A2);

end