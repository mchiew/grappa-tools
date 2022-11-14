function out = apply_weights(input, w)

%   Helper function to apply slice unaliasing grappa weights
%
%   MChiew
%   Nov 2016

%   input is [c,kx,ky,1,t]
%   w is a struct with the kernel and weights ([c,kernel,z])
%   kernel is (kx, ky)

dims    =   size(input);
dims(5) =   size(input,5);
out     =   zeros([dims(1) prod(dims(2:3)) size(w,3) dims(5)]);

%   Pad input boundaries
pad     =   ceil((w.ker-1)/2);
input   =   padarray(input, [0, pad, 0, 0]);

%   Get source and target indices
[src_idx, trg_idx]  =   get_indices(dims(2:3)+2*pad, w.ker);

%   Set periodic boundary condition
input(:, 1:pad(1), :, :, :) =   input(:, end-2*pad(1)+1:end-pad(1), :, :, :);
input(:, end-pad(1)+1:end, :, :, :) =   input(:, pad(1)+1:2*pad(1), :, :, :);
input(:, :, 1:pad(2), :, :) =   input(:, :, end-2*pad(2)+1:end-pad(2), :, :);
input(:, :, end-pad(2)+1:end, :, :) =   input(:, :, pad(2)+1:2*pad(2), :, :);

%   Reshape input to linearize kx, ky
input   =   reshape(input, dims(1), [], 1, dims(5));

for t = 1:dims(5)
    %Select source points for each time point
    src =   reshape(input(:, src_idx, 1, t), [], length(trg_idx));

    %   Apply weights
    for z = 1:w.z
        out(:, :, z, t) =   w.weights(:, :, z)*src;
    end
end

%   Reshape output
out =   reshape(out,[dims(1:3), w.z, dims(5)]);
