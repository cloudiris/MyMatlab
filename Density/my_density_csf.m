function [y_recover x_recover t] = my_density_csf(method, x, y, T, sigma, A, mask, eps, thr, Num_obs, Num_tot, C, sigma_true)

% method: Reconstruction method
% A: Base Matrix(full)
% y: observed signal(full)
% mask: signal entry with non-zero number of observations
% sigma: observed variance
% x_true(x): ground truth sparse vector
% T: sparsity
% thr: throw-away threshold

if (nargin < 8)
    eps = 1;
    thr = 0;
elseif (nargin < 9)
    thr = 0;
elseif (nargin < 10)
    Num_obs = zeros(size(mask), 1);
    Num_tot = zeros(size(mask), 1);
elseif (nargin < 12)
    y_true = zeros(size(y));
end

N = size(x, 1); %sparse vector length
K = size(y, 1); % signal length
x_recover = zeros(N, 1);
y_recover = zeros(K, 1);

% Full;
A_full = A;
y_full = y;
sigma_full = sigma;

% Observed
y = y(mask);
A = A(mask, :);
sigma = sigma(mask);
samplex = 1:size(y, 1);

if (strcmp(method, 'BP') == 1)
% solve by BP
x0 = A'*inv(A*A')*y;
%take epsilon a little bigger than sigma*sqrt(K)!
epsilon = sqrt(sigma' * sigma) * eps;                                                                                                              
tic;
x_BP = l1qc_logbarrier(x0, A, [], y, epsilon, 1);
y_BP = A_full * x_BP;
% [y_BP x_BP] = my_csf(y, samplex', A, T, 'BP', 0); 
t_BP = toc; 
fprintf(1,'BP number of nonzero weights: %d\n',sum(x_BP~=0));
x_recover = x_BP;
y_recover = y_BP;
t = t_BP;

elseif (strcmp(method, 'TBP') == 1)
%solve by throw-away BP(TBP)
tic
if (nargin >= 13) 
    yy = abs(y_true(mask) - y);
    %keep = find(yy(:) <= quantile(yy, thr));
    keep = find(yy(:) <= thr);
else
    %sigma = sigma + ones(size(sigma)) * max(sigma) ./ (sqrt(Num_obs(mask)));
    %keep = find(sigma(:) <= quantile(sigma,thr));
    if (nargin >= 12) 
        sigma = sqrt(sigma.^2  + ones(size(sigma)) * C ./ (Num_obs(mask) .^ 1.5));
    else
        sigma = sqrt(sigma.^2  + ones(size(sigma)) * max(y) .^ 2 ./ (Num_obs(mask) .^ 1.5));
    end
    keep = find(sigma(:) <= thr);
end
x0 = A(keep,:)'*inv(A(keep,:)*A(keep,:)')*y(keep);
epsilon = sqrt(sigma(keep)' * sigma(keep)) * eps;
x_TBP = l1qc_logbarrier(x0, A(keep, :), [], y(keep), epsilon, 1e-3);
y_TBP = A_full * x_TBP;
t_TBP = toc;
fprintf(1,'TBP number of nonzero weights: %d\n',sum(x_TBP~=0));
x_recover = x_TBP;
y_recover = y_TBP;
t = t_TBP;

elseif (strcmp(method, 'DBP') == 1)
%solve by DBP                                                                                                             
tic;
%take epsilon a little bigger than sigma*sqrt(K)! 
if (nargin >= 13) 
    %yy = abs(y_true(mask) - y) + 0.01;
    sigma = sigma_true(mask) ./ sqrt(Num_obs(mask));
else
    %sigma = sigma + ones(size(sigma)) * max(sigma) ./ (sqrt(Num_obs(mask)));
    if (nargin >= 12) 
        sigma = sqrt(sigma.^2  + ones(size(sigma)) * C ./ (Num_obs(mask) .^ 1.5));
    else
        sigma = sqrt(sigma.^2  + ones(size(sigma)) * max(y) .^ 2 ./ (Num_obs(mask) .^ 1.5));
    end
    %keep = find(sigma(:) <= max(sigma));
end
Ab = inv(diag(sigma .^ 2)) * A;
yb = inv(diag(sigma .^ 2)) * y;
epsilon = sqrt(size(y, 1)) * eps;
x0 = Ab'*inv(Ab*Ab')*yb;
x_DBP = l1qc_logbarrier(x0, Ab, [], yb, epsilon, 1e-3);
y_DBP = A_full * x_DBP;
% [y_DBP x_DBP] = my_csf(yb, samplex', Ab, T, 'BP', 0); 
% y_DBP = A * x_DBP;
t_DBP = toc; 
fprintf(1,'DBP number of nonzero weights: %d %f\n',sum(x_DBP~=0), sigma(1));
x_recover = x_DBP;
y_recover = y_DBP;
t = t_DBP;

elseif (strcmp(method, 'TDBP') == 1)
%solve by throw-away DBP(TDBP)
tic
if (nargin >= 13) 
    yy = abs(y_true(mask) - y);
    %keep = find(yy(:) <= quantile(yy, thr));
    sigma = yy + 0.01;
    thr = max(yy);
    keep = find(yy(:) <= thr);
else
    %keep = find(sigma(:) <= quantile(sigma,thr));
    thr = max(sigma);
    keep = find(sigma(:) <= thr);
    sigma = sigma + ones(size(sigma)) * max(sigma) ./ (sqrt(Num_obs(mask)));
    %keep = find(sigma(:) <= max(sigma));
end
Ab = inv(diag(sigma .^ 2)) * A;
yb = inv(diag(sigma .^ 2)) * y;
x0 = Ab(keep,:)'*inv(Ab(keep,:)*Ab(keep,:)')*yb(keep);
epsilon = sqrt(size(yb, 1)) * eps;
x_TDBP = l1qc_logbarrier(x0, Ab(keep, :), [], yb(keep), epsilon, 1e-3);
y_TDBP = A_full * x_TDBP;
t_TDBP = toc;
fprintf(1,'TDBP number of nonzero weights: %d %f\n',sum(x_TDBP~=0), sigma(1));
x_recover = x_TDBP;
y_recover = y_TDBP;
t = t_TDBP;

elseif (strcmp(method, 'LDBP') == 1)
tic
for i = 1:10
    Ab = inv(diag(sigma .^ 2)) * A;
    yb = inv(diag(sigma .^ 2)) * y;
    %keep = find(sigma(:) <= quantile(sigma,thr));
    keep = find(sigma <= thr);
    x0 = Ab(keep,:)'*inv(Ab(keep,:)*Ab(keep,:)')*yb(keep);
    epsilon = sqrt(size(yb, 1)) * eps;
    x_LDBP = l1qc_logbarrier(x0, Ab(keep, :), [], yb(keep), epsilon, 1e-3);
    y_LDBP = A_full * x_LDBP;
    for j = 1:size(y)
        if (sigma(j) > thr && (y(j) - y_LDBP(j) < thr)) 
            sigma(j) = sigma(j) * 0.5 + (y(j) - y_LDBP(j)) * 0.5;
        end
    end
end
t_LDBP = toc;
fprintf(1,'TDBP number of nonzero weights: %d\n',sum(x_LDBP~=0));
x_recover = x_LDBP;
y_recover = y_LDBP;
t = t_LDBP;

elseif (strcmp(method, 'OMP') == 1)
% solve by OMP
tic
[y_OMP x_OMP] = my_csf(y, samplex', A, T, 'OMP', 0);
y_OMP = A_full * x_OMP;
t_OMP = toc;
fprintf(1,'OMP number of nonzero weights: %d\n',sum(x_OMP~=0));
x_recover = x_OMP;
y_recover = y_OMP;
t = t_OMP;

elseif (strcmp(method, 'DOMP') == 1)
% solve by DOMP
tic
if (nargin >= 13) 
    yy = abs(y_true(mask) - y);
    sigma = yy;
else
    %sigma = sigma + ones(size(sigma)) * max(y) ./ (Num_obs(mask) .^ 2);
    if (nargin >= 12) 
        sigma = sqrt(sigma.^2  + ones(size(sigma)) * C ./ (Num_obs(mask) .^ 1.5));
    else
        sigma = sqrt(sigma.^2  + ones(size(sigma)) * max(y) .^ 2 ./ (Num_obs(mask) .^ 1.5));
    end
    %sigma = sqrt(sigma.^2  + ones(size(sigma)) * max(y) .^ 2 ./ (Num_obs(mask) .^ 1.5));
end

[y_DOMP x_DOMP] = my_csf(y, samplex', A, T, 'DOMP', 0, sigma);
y_DOMP = A_full * x_DOMP;
t_DOMP = toc;
fprintf(1,'DOMP number of nonzero weights: %d\n',sum(x_DOMP~=0));
%[bestre co] = my_csf(s, samplex, base, k, method, er)
x_recover = x_DOMP;
y_recover = y_DOMP;
t = t_DOMP;

elseif (strcmp(method, 'TOMP') == 1) 
% solve by TOMP
tic
if (nargin >= 13) 
    yy = abs(y_true(mask) - y);
    %keep = find(yy(:) <= quantile(yy, thr));
    keep = find(yy(:) <= thr);
else
    %keep = find(sigma(:) <= quantile(sigma,thr));
    keep = find(sigma(:) <= thr);
    if (nargin >= 12) 
        sigma = sqrt(sigma.^2  + ones(size(sigma)) * C ./ (Num_obs(mask) .^ 1.5));
    else
        sigma = sqrt(sigma.^2  + ones(size(sigma)) * max(y) .^ 2 ./ (Num_obs(mask) .^ 1.5));
    end
    %sigma = sqrt(sigma.^2  + ones(size(sigma)) * max(y) .^ 2 ./ (Num_obs(mask) .^ 1.5));
    %sigma = sigma + ones(size(sigma)) * max(sigma) ./ (sqrt(Num_obs(mask)));
    %sigma = sigma + ones(size(sigma)) * max(y) ./ (Num_obs(mask) .^ 2);
    %keep = find(sigma(:) <= max(sigma));
end

[y_TOMP x_TOMP] = my_csf(y(keep), samplex(keep)', A, T, 'OMP', 0);
y_TOMP = A_full * x_TOMP;
t_TOMP = toc;
fprintf(1,'TOMP number of nonzero weights: %d\n',sum(x_TOMP~=0));
%[bestre co] = my_csf(s, samplex, base, k, method, er)
x_recover = x_TOMP;
y_recover = y_TOMP;
t = t_TOMP;

elseif (strcmp(method, 'TDOMP') == 1) 
% solve by TDOMP
tic
if (nargin >= 13) 
    yy = abs(y_true(mask) - y);
    %keep = find(yy(:) <= quantile(yy, thr));
    keep = find(yy(:) <= thr);
else
    %keep = find(sigma(:) <= quantile(sigma,thr));
    keep = find(sigma(:) <= thr);
    if (nargin >= 12) 
        sigma = sqrt(sigma.^2  + ones(size(sigma)) * C ./ (Num_obs(mask) .^ 1.5));
    else
        sigma = sqrt(sigma.^2  + ones(size(sigma)) * max(y) .^ 2 ./ (Num_obs(mask) .^ 1.5));
    end
    %sigma = sqrt(sigma.^2  + ones(size(sigma)) * max(y) .^ 2 ./ (Num_obs(mask) .^ 1.5));
    %sigma = sigma + ones(size(sigma)) * max(y) ./ (Num_obs(mask) .^ 2);
    %keep = find(sigma(:) <= max(sigma));
end
[y_TDOMP x_TDOMP] = my_csf(y(keep), samplex(keep)', A, T, 'DOMP', 0, sigma);
y_TDOMP = A_full * x_TDOMP;
t_TDOMP = toc;
fprintf(1,'TDOMP number of nonzero weights: %d\n',sum(x_TDOMP~=0));
%[bestre co] = my_csf(s, samplex, base, k, method, er)
x_recover = x_TDOMP;
y_recover = y_TDOMP;
t = t_TDOMP;

elseif (strcmp(method, 'BCS') == 1)
% solve by BCS
initsigma2 = std(y)^2/1e2;
tic;
[weights,used,sigma2,errbars] = BCS_fast_rvm(A,y,initsigma2,1e-8);
t_BCS = toc;
fprintf(1,'BCS number of nonzero weights: %d\n',length(used));
x_BCS = zeros(N,1); err = zeros(N,1);
x_BCS(used) = weights; err(used) = errbars;
y_BCS = A_full * x_BCS;
x_recover = x_BCS;
y_recover = y_BCS;
t = t_BCS;
end

