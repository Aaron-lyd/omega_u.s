function e = e_NTP(S, T, Z, s, t, z, D, e_tb, e_cb)

% The capital letters of S, T, Z are 3D data and the lower letters of s, t, z are 2d data

lead1 = @(x) reshape(x, [1 size(x)]);

Tppc = ppc_pchip(Z,T); % 3d data on the grid center (g.RC)
Tz = ppc_deriv(Z,Tppc,1); % 3d Tz on the grid center
Tz_2d = ppc_val(Z,Tz,lead1(z)); % interpolate on the 2d surface

Sppc = ppc_pchip(Z,S);
Sz = ppc_deriv(Z,Sppc,1);
Sz_2d = ppc_val(Z,Sz,lead1(z));

Tzz = dfdzz3D(T, Z); % 3d Tzz on the grid center (g.RC)
Tzz_2d = pchipdqn(lead1(z), Z, Tzz); %Tzz on the 2d surface

Szz = dfdzz3D(S, Z);
Szz_2d = pchipdqn(lead1(z), Z, Szz);

% rho_S and rho_z
[rs, rt] = eos_s_t(s, t, z);

sigma_z = -(rs.*Sz_2d + rt.*Tz_2d);

% make the D to be 3D profile to calculate Dz
D_3d = D*ones(size(S));

Dppc = ppc_pchip(Z,D_3d);
Dz_3d = ppc_deriv(Z,Dppc,1);
Dz = ppc_val(Z,Dz_3d,lead1(z));

e = Dz + D*(rt.*Tzz_2d + rs.*Szz_2d)./sigma_z + e_tb + e_cb;
end

