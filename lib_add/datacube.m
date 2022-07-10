function [S, T, P, U, V, R, phi] = datacube(nx,ny,nz,helicity,wall,U_M, periodic, land)

% define some functions
im1 = @(F) circshift(F, [0 +1 0]);
jm1 = @(F) circshift(F, [0 0 +1]);
ip1 = @(F) circshift(F, [0 -1 0]);
jp1 = @(F) circshift(F, [0 0 -1]);

% avergece and backward difference
D_X_f = @(F)  ip1(F) - F;
D_Y_f = @(F)  jp1(F) - F;

D_X_b = @(F)  F - im1(F);
D_Y_b = @(F)  F - jm1(F);

% nx, ny, nz are the size of the datacube
P = (0:100:(nz-1)*100).';
X = linspace(-1,1,nx);
Y = reshape(linspace(-1,1,ny), [1 1 ny]);

latx = 1:nx;
laty = 1:ny;
X_sin = sin(latx*pi/nx);
Y_sin = reshape(sin(laty*pi/ny), [1 1 ny]);

T = repmat((12.5 : -0.5 : 12.5-0.5*(nz-1)).', [1, nx, ny]);
S = repmat((32.1 : .2 : 32.1+0.2*(nz-1)).', [1, nx, ny]);

if periodic
    T = T + X_sin + Y_sin;
    S = S + X_sin + Y_sin;
else
    T = T + X + Y;
    S = S + X + Y;
end

if land
    S(:,wall,:) = nan;
    S(:,:,wall) = nan;
    
    T(:,wall,:) = nan;
    T(:,:,wall) = nan;
end

% Add some helicity
P3 = repmat(P, [1 nx ny]);
R = gsw_rho(S, T, P3);
T = T + helicity * exp(-(X*10).^2 - (Y*10).^2);
S = gsw_SA_from_rho(R, T, P3);
R = gsw_rho(S, T, P3);

%% Build velocity

% streamfunction
phi = (0.5*exp((-X.^2 - Y.^2)/0.1));

% make sure there is no normal flow
% phi(:,1,:) = 0;
% phi(:,end,:) = 0;
% phi(:,:,1) = 0;
% phi(:,:,end) = 0;

if U_M==1
    % create u and v by using backward difference of phi
    U = zeros(nz,nx,ny);
    for i = 1:20
        U(i,:,:) = (1 - 4.5e-2*(i-1))*D_Y_f(phi);
    end
    
    
    V = zeros(nz,nx,ny);
    for i = 1:20
        V(i,:,:) = -(1 - 4.5e-2*(i-1))*D_X_f(phi);
    end
elseif U_M==0
    % create u and v by using backward difference of phi divergence free
    U = zeros(nz,nx,ny);
    for i = 1:20
        U(i,:,:) = (1 - 4.5e-2*(i-1))*D_X_f(phi);
    end
    
    
    V = zeros(nz,nx,ny);
    for i = 1:20
        V(i,:,:) = (1 - 4.5e-2*(i-1))*D_Y_f(phi);
    end
    
elseif U_M==2
    
    U = zeros(nz,nx,ny);
    for i = 1:20
        U(i,:,:) = (100*sin(i/2)+0.6)*(phi - circshift(phi, [0, 0 +1]));
    end
    
    
    V = zeros(nz,nx,ny);
    for i = 1:20
        V(i,:,:) = -(100*cos(i/2)+0.6)*(phi - circshift(phi, [0, +1 0]));
    end
    
elseif U_M==3
    % artificail velocity
    
    U = ones(nz,nx,ny,'double');
    V = zeros(nz,nx,ny,'double');
    
elseif U_M==4
    
    U = rand(nz,nx,ny,'double');
    V = rand(nz,nx,ny,'double');
    
elseif U_M==5
    
    U = ones(nz,nx,ny);
    V = ones(nz,nx,ny);
    
end

end