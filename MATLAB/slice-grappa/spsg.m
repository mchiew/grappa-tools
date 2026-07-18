function out = spsg(input, calib, kernel)

%   Split-slice GRAPPA Multi-band Recon
%   Based on Cauley et al., MRM 2014
%
%   MChiew
%   Nov 2016

%   input is [c,kx,ky,1,t]
%   calib is [c,kx,ky,z]
%   kernel is (kx, ky)

w   =   weights_spsg(calib, kernel);
out =   apply_weights(input, w);
