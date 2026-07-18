function [src, trg] = get_indices(dims, kernel)

%   Helper function to compute source, target and kernel indices
%   
%   MChiew
%   Nov 2016

%   dims is (nx, ny)
%   kernel is (kx, ky), where kx and ky should be odd

pad =   ceil((kernel-1)/2);
ks  =   prod(kernel);

%   Find calib boundary padding
kx  =   1+pad(1):dims(1)-pad(1);
ky  =   1+pad(2):dims(2)-pad(2);

%   Find relative kernel indices
mask    =   false(dims);
mask(1:kernel(1), 1:kernel(2))  =   true;
k_idx   =   find(mask);
k_idx   =   k_idx - k_idx((ks+1)/2);

%   Find target linear indices
mask    =   false(dims);
mask(kx, ky)    =   true;
trg     =   find(mask);

%   Find source linear indices
src =   bsxfun(@plus, k_idx, trg');
