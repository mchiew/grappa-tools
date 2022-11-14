%   grappa.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           data    -   (nc, nx, ny, nz, m]) complex undersampled k-space data
%                       will also loop across extra dimension m
%           calib   -   (nc, cx, cy, cz) complex calibration k-space data
%           R       -   [Rx, Ry] or [Rx, Ry, Rz] acceleration factors
%           kernel  -   [kx, ky] or [kx, ky, kz] kernel size 
%           tol     -   singular value cutoff threshold for kernel weight
%                       training, relative to s(1), defaults to pinv default
%
%   output:
%           recon   -   (nc, nx, ny, nz) complex reconstructed k-space data

function data = grappa(data, calib, R, kernel, tol)

%% Use default pinv tolerance if not supplied
if nargin < 5
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

%%  Prepare masks and zero-pad data
pad     =   floor(R.*kernel/2);
mask    =   padarray(data~=0, [0 pad]);
data    =   padarray(data,    [0 pad]);
loop    =   size(data,5);
offset  =   numel(data(:,:,:,:,1));

%%  Loop over all possible kernel types
for type = 1:prod(R(2:end))-1

    %   Collect source and target calibration points for weight estimation   
    [src, trg]  =   grappa_get_indices(kernel, true(size(calib)), pad, R, type);

    %   Perform weight estimation    
    weights     =   calib(trg)*pinv_reg(calib(src));
    
    %   Loop over extra dimension if they exist
    for m = 1:loop               
        
        %   Collect source points in under-sampled data for weight application    
        [src, trg]  =   grappa_get_indices(kernel, mask(:,:,:,:,m), pad, R, type, (m-1)*offset);
        
        %   Apply weights to reconstruct missing data    
        data(trg)   =   weights*data(src);         
    
    end
    
end

%%  Un-pad reconstruction to get original image size back
data   =   data(:,pad(1)+1:size(data,2)-pad(1), pad(2)+1:size(data,3)-pad(2), pad(3)+1:size(data,4)-pad(3),:);
