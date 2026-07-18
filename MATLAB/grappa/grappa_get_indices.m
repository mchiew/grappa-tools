%   grappa_get_indices.m
%   mchiew@fmrib.ox.ac.uk
%
%   inputs: 
%           kernel  -   [sx, sy, sz] kernel size in each dimension
%           samp    -   (c, nx, ny, nz) sampling mask
%           pad     -   [pad_x, pad_y, pad_z] size of padding in each direction 
%           type    -   (scalar, must be < R) indicates which of the R(2)*R(3)-1 kernels
%                       you are trying to index over
%           offset  -   additional index offset that gets added to src,trg
%
%   output:
%           src     -   linear indices for all source points (c*sx*sy*sz, all possible targets)
%           trg     -   linear indices for all the target points (c, all possible targets)

function [src, trg] = grappa_get_indices(kernel, samp, pad, R, type, offset)

%   Offset is optional, 0 by default
if nargin < 6
    offset  =   0;
end

%   Get dimensions
[nc,dx,dy,dz]    =   size(samp);

%   Make sure the under-sampling is in y and z only
%   There are a few things here that require that assumption
if R(1) > 1
    error('x-direction must be fully sampled');
end

%   Make sure the type parameter makes sense
%   It should be between 1 and R(2)*R(3)-1 (inclusive)
if type > prod(R(2:3))-1
    error('Type parameter is inconsistent with R');
end

%   Find the limits of all possible target points given padding
kx  =   1+pad(1):dx-pad(1);
ky  =   1+pad(2):dy-pad(2);
kz  =   1+pad(3):dz-pad(3);

%%  Compute indices for a single coil

%   Find relative indices for kernel source points
mask    =   false(dx,dy,dz);
mask(1:R(1):R(1)*kernel(1), 1:R(2):R(2)*kernel(2), 1:R(3):R(3)*kernel(3))    =   true;
k_idx   =   reshape(find(mask),[],1);

%   Find the index for the desired target point (depends on type parameter)
mask    =   false(dx,dy,dz);
[yy,zz] =   ind2sub(R(2:3),type+1);
mask(R(1)*ceil(kernel(1)/2), R(2)*(ceil(kernel(2)/2)-1)+yy, R(3)*(ceil(kernel(3)/2)-1)+zz)    =   true;
k_trg   =   reshape(find(mask),[],1);

%   Subtract the target index from source indices
%   to get relative linear indices for all source points
%   relative to the target point (index 0, target position)
k_idx   =   k_idx - k_trg;

%   Find all possible target indices
mask    =   false(dx,dy,dz);
mask(kx,ky,kz) =   squeeze(circshift(samp(1,kx,ky,kz),[0 0 yy-1 zz-1]));
trg     =   reshape(find(mask),1,[]);

%   Find all source indices associated with the target points in trg
src =   bsxfun(@plus, k_idx, trg);

%%  Now replicate indexing over all coils

%   Final shape of trg should be (#coils, all possible target points)
trg =   bsxfun(@plus, (trg-1)*nc+1, (0:nc-1)') + offset;

%   Final shape of src should be (#coils*sx*sy, all possible target points)
src =   bsxfun(@plus, (src(:)'-1)*nc+1, (0:nc-1)');
src =   reshape(src,[], size(trg,2)) + offset;
