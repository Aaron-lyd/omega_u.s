function [x_venm_x,x_venm_y] = x_VENM(S, T, X, s ,t ,x, OPTS)

[nk,ni,nj] = size(S);
% code generation
ntp_bottle_to_cast_codegen(nk,OPTS);

TOL_DENS = 1e-7;
drdp = 4e-03; % Approximate derivative of in-situ density w.r.t. pressure [kg m^-3 dbar^-1]
X_TOL = TOL_DENS / (drdp * 2); % tolerance in pressure [dbar] during vertical solve.  Factor of 2 for good measure
 
SppX = ppc_linterp(X, S);
TppX = ppc_linterp(X, T);

x_venm_x = zeros(size(x));
x_venm_y = zeros(size(x));

BotK = squeeze(sum(isfinite(S), 1));
for j = 2:nj
    jm1 = j-1;
    for i = 1:ni
        k = BotK(i,j);
        im1 = mod(i-2, ni) + 1;
        x_venm_x(i,j) = ntp_bottle_to_cast_mex(SppX(:,:,im1,j), TppX(:,:,im1,j), X, k, s(i,j), t(i,j), x(i,j), X_TOL);

        x_venm_y(i,j) = ntp_bottle_to_cast_mex(SppX(:,:,i,jm1), TppX(:,:,i,jm1), X, k, s(i,j), t(i,j), x(i,j), X_TOL);
        
    end
end