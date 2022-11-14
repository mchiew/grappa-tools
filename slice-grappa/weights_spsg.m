function w = weights_spsg(calib, kernel, lambda)

%   Helper function to compute split slice grappa kernel weights
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

for z = 1:dims(4)
    %   Select target points for each slice
    trg =   calib(:, trg_idx, z);

    %   Select source points for target slice
    src =   reshape(calib(:, src_idx, z), [], length(trg_idx));

    %   Select all non-target slice locations
    zz  =   setdiff(1:dims(4), z);

    %   Generate zeroed target points
    trg2=   repmat(zeros(size(trg)), 1, length(zz));

    %   Select source points for non-target slice
    src2=   reshape(calib(:, src_idx, zz), [], length(trg_idx)*length(zz));

    %   Get (regularised) pseudoinverse
    m   =   [src src2]'*pinv([src src2]*[src src2]' + norm(src)*lambda*eye(size(src,1)));

    %   Fit weights
    w.weights(:, :, z)  =   [trg trg2]*m;
end
