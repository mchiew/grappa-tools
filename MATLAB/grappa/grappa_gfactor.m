%   grappa_gfactor.m
%   mark.chiew@ndcn.ox.ac.uk
%
%   inputs: 
%           data    -   (nc, nx, ny, nz]) complex undersampled k-space data
%                       recon will be applied to any extra dimensions
%           calib   -   (nc, cx, cy, cz) complex calibration k-space data
%           noise   -   (c,c) coil noise covariance matrix
%           R       -   [Rx, Ry] or [Rx, Ry, Rz] acceleration factors
%           kernel  -   [kx, ky] or [kx, ky, kz] kernel size 
%           tol     -   singular value cutoff threshold for kernel weight
%                       training, relative to s(1), defaults to pinv
%
%   output:
%           data    -   image-domain unaliased reconstructed dataset
%           g       -   (nx ny, nz) sos-combined g-factor

function [data, g]  =   grappa_gfactor(data, calib, noise, R, kernel, tol)

%% Use default pinv tolerance if not supplied
if nargin < 6
    pinv_reg = @pinv;
else
    pinv_reg = @(A)pinv(A, tol*norm(A,2));
end

%% Determine whether this is a 1D or 2D GRAPPA problem
if numel(R) == 2
    R(3)        =   1;
end
if numel(kernel) == 2
    kernel(3)   =   1;
end

%% Calib dimensions
Nc  =   size(calib,1);
Nx  =   size(calib,2);
Ny  =   size(calib,3);
Nz  =   size(calib,4);

%% Output dimensions
Mx  =   size(data,2);
My  =   size(data,3);
Mz  =   size(data,4);

%% Setup 
W   =   zeros([Nc, Nx*Ny*Nz, Nc]);

mask    =   false([Nc Nx Ny Nz]); 
mask(:,floor(end/2+1),floor(end/2+1),floor(end/2+1)) =   true;

%%  Loop over all possible kernel types
for type = 1:prod(R(2:end))-1

    %   Collect source and target calibration points for weight estimation   
    [src, trg]  =   grappa_get_indices(kernel, true([Nc Nx Ny Nz]), floor(R.*kernel/2), R, type);

    %   Perform weight estimation    
    weights     =   calib(trg)*pinv_reg(calib(src));
    
    %   Get indices for building kernel image
    [yy,zz]     =   ind2sub(R(2:3),type+1);
    [src, ~]    =   grappa_get_indices(kernel, circshift(mask, [0 0 1-yy 1-zz]), floor(R.*kernel/2), R, type);
    
    %   Reshape data
    weights     =   reshape(weights, Nc, Nc, []);
    
    %   Populate kernel image k-space
    for c = 1:Nc
        idx     =   src+(c-1)*prod([Nc Nx Ny Nz]);
        W(idx)  =   weights(c,:,:);
    end
end

%% Flip, pad and inverse Fourier Transform kernels
W   =   reshape(W, [Nc, Nx, Ny, Nz, Nc]);
W   =   circshift(flip(W,2),1,2);
W   =   circshift(flip(W,3),1,3);
W   =   circshift(flip(W,4),1,4);
W   =   padarray(W, [0 ceil(([Mx,My,Mz]-[Nx,Ny,Nz])/2) 0], 'pre');
W   =   padarray(W, [0 floor(([Mx,My,Mz]-[Nx,Ny,Nz])/2) 0], 'post');
W   =   permute(W, [5,1,2,3,4]);
W(:,:,floor(end/2+1),floor(end/2+1),floor(end/2+1)) = eye(Nc);
W   =   ifftdim(W,3:5)*sqrt(Mx*My*Mz);

%% Apply image-kernels to unalias data
data    =   shiftdim(ifftdim(data,2:4),-1);;
data    =   squeeze(sum(W.*data,2));

%% Compute g-factor
W   =   reshape(W, Nc, Nc, []);
data=   reshape(conj(data), Nc, 1, []);
W   =   squeeze(sum(data.*W,1));
W   =   sum(W.*(conj(noise*W)),1);
p   =   sum(squeeze(data).*conj(noise*squeeze(data)),1);
g   =   reshape(sqrt(abs(W)./abs(p)), [Mx,My,Mz])/prod(R);
data=   reshape(conj(data), Nc, Mx, My, Mz);
