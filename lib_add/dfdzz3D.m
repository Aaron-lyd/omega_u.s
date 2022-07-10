function dfdzz = dfdzz3D(f, Z)
    % Given f = f(z, x, y) where z takes values given by Z, evaluate d2f/fz2
    % f = a 3D field of data.
    % Z = z coordinate for first dimension of f

    assert(size(f, 1) == length(Z), 'Z does not match size of 1st dim of f');

    Z = Z(:);
    Zm1 = [NaN; Z(1:end-1)];
    Zp1 = [Z(2:end); NaN];

    fm1 = circshift(f, [+1, 0, 0]);
    fp1 = circshift(f, [-1, 0, 0]);

    % Use Lagrange Interpolating Polynomial.
    dfdzz =         fm1 .* 2 ./ ((Zm1 - Z  ) .* (Zm1 - Zp1)) ;
    dfdzz = dfdzz + f   .* 2 ./ ((Z - Zm1  ) .* (Z - Zp1  )) ;
    dfdzz = dfdzz + fp1 .* 2 ./ ((Zp1 - Zm1) .* (Zp1 - Z  )) ;
end