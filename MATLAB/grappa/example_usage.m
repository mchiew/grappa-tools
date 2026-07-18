%% load sensitivities
load('sens.mat')

%% generate single slice test image
x = phantom(64);

%% multiply through sensitivity maps
x = x.*sens;

%% Fourier transform input
y = fftshift(fft(fftshift(fft(x,[],1),1),[],2),2);

%% generate noisy slice aliased input, and calibration data
%  need to permute because grappa expect coils as first dim
calib = permute(y,[4,1,2,3]);
calib = calib(:,21:44,21:44);

input = permute(y,[4,1,2,3]);
noisy = input + (randn(size(input)) + 1j*randn(size(input)))/10;
noisy(:,:,1:2:end) = 0;

%% perform grappa reconstruction
out   = grappa(noisy, calib, [1,2], [3,4]);

%% inverse Fourier transform and sum-of-squares combine to show outputs
img_grappa = squeeze(sqrt(sum(abs(ifft(ifftshift(ifft(ifftshift(out,2),[],2),3),[],3)),1)));
img_aliased = squeeze(sqrt(sum(abs(ifft(ifftshift(ifft(ifftshift(noisy,2),[],2),3),[],3)),1)));
