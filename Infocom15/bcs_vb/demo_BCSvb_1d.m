% Demo for BCSvb on 1-dimensional signal
% Lihan He, Mar. 3, 2009

%----------------
% Generate signal
%----------------

% Generate original step signal -- signal0
% signal0: random generated step signal with random step position, step
% width and step amplitude
Nsignal=512;    % signal length
Nstep=10;       % number of steps in the entire signal 
maxWstep=30;    % max width of each step 
Astep=1;        % approximate amplitude region of each step 
q = randperm(Nsignal);
Tstep=q(1:Nstep);       % random step position
Wstep=round(rand(1,Nstep)*maxWstep);    % random step width
Wstep_l=Tstep-round(Wstep/2);
Wstep_r=Wstep_l+Wstep;
signal0=zeros(1,Nsignal);
for i=1:Nstep
    idx=max(1,Wstep_l(i)):min(Nsignal,Wstep_r(i));
    signal0(idx)=signal0(idx)+Astep*randn*ones(1,length(idx));
end

% 1D wavelet transformation and sparse signal -- theta0 
DecLevel=6;
WaveletName='db1';
[C0, S0] = wavedec(signal0, DecLevel, WaveletName);
theta0=C0(:);

%----------------------------------
% Projection matrix and observation
%----------------------------------

% Generate projection matrix Phi, N x M matrix
N = 150;    % number of CS measurements
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

signal=waverec(theta', S0, WaveletName);
ERR=sqrt(sum((signal-signal0).^2,1)/sum(signal0.^2,1));

% ----
% plot
% ----

figure, subplot(2,2,1), plot(signal0,'r'); 
axis([1,M,1.1*min(signal0), 1.1*max(signal0)]); title('Original Signal')
subplot(2,2,2); plot(theta0,'r')
axis([1,M,1.1*min(theta0), 1.1*max(theta0)]); title('Original Sparse Coefficients')
subplot(2,2,3); plot(theta);
axis([1,M,1.1*min(theta0), 1.1*max(theta0)]); title('Reconstructed Sparse Coefficients')
subplot(2,2,4); plot(signal);
axis([1,M,1.1*min(signal0), 1.1*max(signal0)]); title(['Recovered Signal, rela err = ',num2str(ERR)])
