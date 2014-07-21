
N = 512; % sparse vector length
T = 5;  % sparsity
K = 400; % signal length

[x y_true y sigma A mask] = signal_generator(N, T, K);

A_true = A;
y_obs = y;
y = y(mask);
A = A(mask, :);
sigma = sigma(mask);
samplex = 1:size(y, 1);

% solve by BP
x0 = A'*inv(A*A')*y;
%take epsilon a little bigger than sigma*sqrt(K)!
epsilon =  sqrt(sigma' * sigma);                                                                                                              
tic;
x_BP = l1qc_logbarrier(x0, A, [], y, epsilon, 1);
y_BP = A_true * x_BP;
% [y_BP x_BP] = my_csf(y, samplex', A, T, 'BP', 0); 
t_BP = toc; 
fprintf(1,'BP number of nonzero weights: %d\n',sum(x_BP~=0));

%solve by throw-away BP(TBP)
tic
keep = find(sigma(:) <= quantile(sigma,0.3));
x0 = A(keep,:)'*inv(A(keep,:)*A(keep,:)')*y(keep);
epsilon = sqrt(sigma(keep)' * sigma(keep));
x_TBP = l1qc_logbarrier(x0, A(keep, :), [], y(keep), epsilon, 1e-3);
y_TBP = A_true * x_TBP;
t_TBP = toc;
fprintf(1,'TBP number of nonzero weights: %d\n',sum(x_TBP~=0));

%solve by DBP                                                                                                             
tic;
%take epsilon a little bigger than sigma*sqrt(K)!   
Ab = inv(diag(sigma .^ 2)) * A;
yb = inv(diag(sigma .^ 2)) * y;
epsilon = sqrt(size(y, 1)) * 5;
x0 = Ab'*inv(Ab*Ab')*yb;
x_DBP = l1qc_logbarrier(x0, Ab, [], yb, epsilon, 1e-3);
y_DBP = A_true * x_DBP;
% [y_DBP x_DBP] = my_csf(yb, samplex', Ab, T, 'BP', 0); 
% y_DBP = A * x_DBP;
t_DBP = toc; 
fprintf(1,'DBP number of nonzero weights: %d\n',sum(x_DBP~=0));

%solve by throw-away DBP(TDBP)
tic
x0 = Ab(keep,:)'*inv(Ab(keep,:)*Ab(keep,:)')*yb(keep);
epsilon = sqrt(size(yb, 1)) * 5;
x_TDBP = l1qc_logbarrier(x0, Ab(keep, :), [], yb(keep), epsilon, 1e-3);
y_TDBP = A_true * x_TDBP;
t_TDBP = toc;
fprintf(1,'TDBP number of nonzero weights: %d\n',sum(x_TDBP~=0));

% solve by OMP
tic
[y_OMP x_OMP] = my_csf(y, samplex', A, T, 'OMP', 0);
y_OMP = A_true * x_OMP;
t_OMP = toc;
fprintf(1,'OMP number of nonzero weights: %d\n',sum(x_OMP~=0));

% solve by DOMP
tic
[y_DOMP x_DOMP] = my_csf(y, samplex', A, T, 'DOMP', 0, sigma);
y_DOMP = diag(sigma) * A_true * x_DOMP;
t_DOMP = toc;
fprintf(1,'DOMP number of nonzero weights: %d\n',sum(x_DOMP~=0));
%[bestre co] = my_csf(s, samplex, base, k, method, er)

% solve by TDOMP
tic
[y_TDOMP x_TDOMP] = my_csf(y(keep), samplex(keep)', A, T, 'DOMP', 0, sigma);
y_TDOMP = diag(sigma) * A_true * x_TDOMP;
t_TDOMP = toc;
fprintf(1,'TDOMP number of nonzero weights: %d\n',sum(x_TDOMP~=0));
%[bestre co] = my_csf(s, samplex, base, k, method, er)

% solve by BCS
initsigma2 = std(y)^2/1e2;
tic;
[weights,used,sigma2,errbars] = BCS_fast_rvm(A,y,initsigma2,1e-8);
t_BCS = toc;
fprintf(1,'BCS number of nonzero weights: %d\n',length(used));
x_BCS = zeros(N,1); err = zeros(N,1);
x_BCS(used) = weights; err(used) = errbars;
y_BCS = A_true * x_BCS;

% reconstruction error
E_BP = norm(x-x_BP)/norm(x);
E_BCS = norm(x-x_BCS)/norm(x);
E_DBP = norm(x-x_DBP)/norm(x);
E_TDBP = norm(x-x_TDBP)/norm(x);
E_TBP = norm(x-x_TBP)/norm(x);

E_OMP = norm(x-x_OMP)/norm(x);
E_DOMP = norm(x-x_DOMP)/norm(x);
E_TDOMP = norm(x-x_TDOMP)/norm(x);

Ey_BP = norm(y_true-y_BP)/norm(y_true);
Ey_BCS = norm(y_true-y_BCS)/norm(y_true);
Ey_DBP = norm(y_true-y_DBP)/norm(y_true);
Ey_TDBP = norm(y_true-y_TDBP)/norm(y_true);
Ey_TBP = norm(y_true-y_TBP)/norm(y_true);

Ey_OMP = norm(y_true-y_OMP)/norm(y_true);
Ey_DOMP = norm(y_true-y_DOMP)/norm(y_true);
Ey_TDOMP = norm(y_true-y_TDOMP)/norm(y_true);


% figure
% subplot(3,1,1); plot(x); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(a) Original Signal']);
% subplot(3,1,2); plot(x_BP); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(b) Reconstruction with BP, K=' num2str(K)]);
% subplot(3,1,3); errorbar(x_BCS,err); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(c) Reconstruction with BCS, K=' num2str(K)]); box on;

figure;
subplot(4,1,1);
plot(y_true);
hold on;
plot(y_BP, 'r');
title(['(a) BP']);

subplot(4,1,2);
plot(y_true);
hold on;
plot(y_TBP, 'g');
title(['(b) Throw-away BP']);


subplot(4,1,3);
plot(y_true);
hold on;
plot(y_DBP, 'm');
title(['(c) Density BP']);

subplot(4,1,4);
plot(y_true);
hold on;
plot(y_TDBP, 'm');
title(['(d) Throw-away Density BP']);

figure;
subplot(4,1,1);
plot(y_true);
hold on;
plot(y_OMP, 'r');
title(['(a) OMP']);

subplot(4,1,2);
plot(y_true);
hold on;
plot(y_DOMP, 'g');
title(['(b) density OMP']);


subplot(4,1,3);
plot(y_true);
hold on;
plot(y_TDOMP, 'm');
title(['(c) Throw-away Density OMP']);

subplot(4,1,4);
plot(y_true);
hold on;
plot(y_BCS, 'm');
title(['(c) Bayeian Compressive Sensing']);ff

figure;
subplot(4,1,1);
plot(x);
hold on;
plot(x_BP, 'r');
title(['(a) BP']);

subplot(4,1,2);
plot(x);
hold on;
plot(x_TBP,'m');
title(['(b) Throw-away BP']);

subplot(4,1,3);
plot(x);
hold on;
plot(x_DBP, 'g');
title(['(c) Density BP']);

subplot(4,1,4);
plot(x);
hold on;
plot(x_TDBP,'r');
title(['(d) Throw-away Density BP']);

figure;
subplot(4,1,1);
plot(x);
hold on;
plot(x_OMP, 'r');
title(['(a) OMP']);

subplot(4,1,2);
plot(x);
hold on;
plot(x_DOMP,'m');
title(['(b) Density OMP']);

subplot(4,1,3);
plot(x);
hold on;
plot(x_TDOMP, 'g');
title(['(c) Throw-away Density OMP']);

subplot(4,1,4);
plot(x);
hold on;
plot(x_BCS,'r');
title(['(d) Bayesian Compressive Sensing']);

disp(['BP: ||I_hat-I||/||I|| = ' num2str(E_BP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_BP) ', time = ' num2str(t_BP) ' secs']);
disp(['TBP: ||I_hat-I||/||I|| = ' num2str(E_TBP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_TBP) ', time = ' num2str(t_TBP) ' secs']);
disp(['DBP: ||I_hat-I||/||I|| = ' num2str(E_DBP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_DBP) ', time = ' num2str(t_DBP) ' secs']);
disp(['TDBP: ||I_hat-I||/||I|| = ' num2str(E_TDBP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_TDBP) ', time = ' num2str(t_TDBP) ' secs']);

disp(['OMP: ||I_hat-I||/||I|| = ' num2str(E_OMP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_OMP) ', time = ' num2str(t_OMP) ' secs']);
disp(['DOMP: ||I_hat-I||/||I|| = ' num2str(E_DOMP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_DOMP) ', time = ' num2str(t_DOMP) ' secs']);
disp(['TDOMP: ||I_hat-I||/||I|| = ' num2str(E_TDOMP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_TDOMP) ', time = ' num2str(t_TDOMP) ' secs']);
disp(['BCS: ||I_hat-I||/||I|| = ' num2str(E_BCS) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_BCS) ', time = ' num2str(t_BCS) ' secs']);


