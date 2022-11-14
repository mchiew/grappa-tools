function w = weights_sg(calib, kernel, lambda)

%   Helper function to compute slice grappa kernel weights
%   
%   MChiew
%   Nov 2016

%   calib is [c,kx,ky,z]
%   kernel is (kx, ky), where kx and ky should be odd
%   lambda is an optional regularisation parameter relative to norm(src)

if nargin < 3
    lambda = 0;
end

dims    =   size(calib);
w.ker   =   kernel;
w.z     =   dims(4);

%   Get source and target indices
[src_idx, trg_idx]  =   get_indices(dims(2:3), w.ker);

%   Reshape calib to linearize kx, ky
calib   =   reshape(calib, dims(1), [], dims(4));

%   Select source points on summed calibration data
src =   reshape(sum(calib(:, src_idx, :),3), [], length(trg_idx));

for z = 1:dims(4)
    %   Select target points for each slice
    trg =   calib(:, trg_idx, z);
    
    %   Get (regularised) pseudoinverse
    m   =   src'*pinv(src*src' + norm(src)*lambda*eye(size(src,1)));

    %   Fit weights
    w.weights(:, :, z)  =   trg*m;
end
