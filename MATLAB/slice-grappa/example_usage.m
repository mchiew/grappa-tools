%% load sensitivities
load('sens.mat')

%% generate 4-slice test image
x        = zeros(64,64,4);
x(:,:,1) = phantom(64);
x(:,:,2) = permute(x(:,:,1),[2,1,3]);
x(:,:,3) = flip(x(:,:,1),1);
x(:,:,4) = flip(x(:,:,2),2);

%% multiply through sensitivity maps
x = x.*sens;

%% Fourier transform input
y = fft(fft(x,[],1),[],2);

%% generate noisy slice aliased input, and calibration data
%  need to permute because sg/spsg expect coils as first dim
calib = permute(y,[4,1,2,3]);
input = permute(sum(y,3),[4,1,2,3]);
noisy = input + (randn(size(input)) + 1j*randn(size(input)))/1000;

%% perform slice grappa (sg) or split-slice grappa (spsg)
out_sg = sg(noisy, calib, [5,5]);
out_spsg = spsg(noisy, calib, [5,5]);

%% inverse Fourier transform
out_sg = ifft(ifft(out_sg,[],2),[],3);
out_spsg = ifft(ifft(out_spsg,[],2),[],3);
