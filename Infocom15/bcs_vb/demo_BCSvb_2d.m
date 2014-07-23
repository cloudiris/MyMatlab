% Demo for BCSvb on 2-dimensional 128 x 128 signal
% Source for the original image "cameraman": http://decsai.ugr.es/cvg/dbimagenes/
% Lihan He, Mar. 3, 2009

%----------------
% Generate signal
%----------------

% Load original image -- image0
% image0: original 2-d image, cut from "cameraman" with rows 21:148 and columns 71:198
I=imread('cameraman.pgm');
I=double(I);
image0=I(21:148,71:198);

% 2D wavelet transformation and sparse signal -- theta0 
DecLevel=4;
WaveletName='db1';
[C0, S0] = wavedec2(image0, DecLevel, WaveletName);
theta0=C0(:);

%----------------------------------
% Projection matrix and observation
%----------------------------------

% Generate projection matrix Phi, N x M matrix
N = 5000;    % number of CS measurements
M = length(theta0);
Phi = randn(N,M);
Phi = Phi./repmat(sqrt(sum(Phi.^2,1)),[N,1]);	

% CS observations
v = Phi*theta0;

%----------------------
% CS inversion by BCSvb
%----------------------

theta = BCSvb(Phi, v);
% "[]" can be used as input arguments for default values, e.g.,
% theta = BCSvb(Phi, v, [], [], [], 1);

%---------------------
% Reconstruction error
% --------------------

image=waverec2(theta', S0, WaveletName);
ERR=sqrt(sum(sum((image-image0).^2,1),2)/sum(sum(image0.^2,1),2));

% ----
% plot
% ----

figure, subplot(2,2,1), imagesc(image0); colormap(gray); 
axis square; title('Original Image')
subplot(2,2,2); plot(theta0,'r')
axis([1,M,1.1*min(theta0), 1.1*max(theta0)]); title('Original Sparse Coefficients')
subplot(2,2,3); plot(theta);
axis([1,M,1.1*min(theta0), 1.1*max(theta0)]); title('Reconstructed Sparse Coefficients')
subplot(2,2,4); imagesc(image); colormap(gray);
axis square; title(['Recovered Image, rela err = ',num2str(ERR)])
