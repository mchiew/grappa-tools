function out = sg(input, calib, kernel)

%   Slice-GRAPPA

%   MChiew
%   Nov 2016

%   input is [c,kx,ky,1,t]
%   calib is [c,kx,ky,z]
%   kernel is (kx, ky)

w   =   weights_sg(calib, kernel);
out =   apply_weights(input, w);
