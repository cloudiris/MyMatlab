
N = 512; % sparse vector length
T = 5;  % sparsity
K = 400; % signal length

[x y_true y sigma A mask Num_obs Num_tot] = signal_readin_traffic(0.3);

N = size(x, 1);
T = 10;
K = size(y, 1);

[y_BP x_BP t_BP] = my_density_csf('BP', x, y, T, sigma, A, mask, 0.2);
[y_DBP x_DBP t_DBP] = my_density_csf('DBP', x, y, T, sigma, A, mask, 0.02, 15, Num_obs, Num_tot);
[y_TBP x_TBP t_TBP] = my_density_csf('TBP', x, y, T, sigma, A, mask, 0.02, 15, Num_obs, Num_tot);
[y_TDBP x_TDBP t_TDBP] = my_density_csf('TDBP', x, y, T, sigma, A, mask, 0.015, 20, Num_obs, Num_tot, y_true);
[y_OMP x_OMP t_OMP] = my_density_csf('OMP', x, y, T, sigma, A, mask);
[y_TOMP x_TOMP t_TOMP] = my_density_csf('TOMP', x, y, T, sigma, A, mask, 0, 10, Num_obs, Num_tot);
[y_DOMP x_DOMP t_DOMP] = my_density_csf('DOMP', x, y, T, sigma, A, mask, 0, 10, Num_obs, Num_tot);
[y_TDOMP x_TDOMP t_TDOMP] = my_density_csf('TDOMP', x, y, T, sigma, A, mask, 0, 10, Num_obs, Num_tot);
[y_BCS x_BCS t_BCS] = my_density_csf('BCS', x, y, T, sigma, A, mask);

% reconstruction error
E_BP = norm(x-x_BP)/norm(x);
E_BCS = norm(x-x_BCS)/norm(x);
E_DBP = norm(x-x_DBP)/norm(x);
E_TDBP = norm(x-x_TDBP)/norm(x);
E_TBP = norm(x-x_TBP)/norm(x);

E_OMP = norm(x-x_OMP)/norm(x);
E_DOMP = norm(x-x_DOMP)/norm(x);
E_TOMP = norm(x-x_TOMP)/norm(x);
E_TDOMP = norm(x-x_TDOMP)/norm(x);

Ey_BP = norm(y_true-y_BP)/norm(y_true);
Ey_BCS = norm(y_true-y_BCS)/norm(y_true);
Ey_DBP = norm(y_true-y_DBP)/norm(y_true);
Ey_TDBP = norm(y_true-y_TDBP)/norm(y_true);
Ey_TBP = norm(y_true-y_TBP)/norm(y_true);

Ey_OMP = norm(y_true-y_OMP)/norm(y_true);
Ey_DOMP = norm(y_true-y_DOMP)/norm(y_true);
Ey_TOMP = norm(y_true-y_TOMP)/norm(y_true);
Ey_TDOMP = norm(y_true-y_TDOMP)/norm(y_true);


% figure
% subplot(3,1,1); plot(x); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(a) Original Signal']);
% subplot(3,1,2); plot(x_BP); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(b) Reconstruction with BP, K=' num2str(K)]);
% subplot(3,1,3); errorbar(x_BCS,err); axis([1 N -max(abs(x))-0.2 max(abs(x))+0.2]); title(['(c) Reconstruction with BCS, K=' num2str(K)]); box on;

% figure;
% subplot(4,1,1);
% plot(y_true);
% hold on;
% plot(y_BP, 'r');
% title(['(a) BP']);
% 
% subplot(4,1,2);
% plot(y_true);
% hold on;
% plot(y_TBP, 'g');
% title(['(b) Throw-away BP']);
% 
% 
% subplot(4,1,3);
% plot(y_true);
% hold on;
% plot(y_DBP, 'm');
% title(['(c) Density BP']);
% 
% subplot(4,1,4);
% plot(y_true);
% hold on;
% plot(y_TDBP, 'm');
% title(['(d) Throw-away Density BP']);
% 
% figure;
% subplot(4,1,1);
% plot(y_true);
% hold on;
% plot(y_OMP, 'r');
% title(['(a) OMP']);
% 
% subplot(4,1,2);
% plot(y_true);
% hold on;
% plot(y_DOMP, 'g');
% title(['(b) density OMP']);
% 
% 
% subplot(4,1,3);
% plot(y_true);
% hold on;
% plot(y_TDOMP, 'm');
% title(['(c) Throw-away Density OMP']);
% 
% subplot(4,1,4);
% plot(y_true);
% hold on;
% plot(y_BCS, 'm');
% title(['(c) Bayeian Compressive Sensing']);
% 
% figure;
% subplot(4,1,1);
% plot(x);
% hold on;
% plot(x_BP, 'r');
% title(['(a) BP']);
% 
% subplot(4,1,2);
% plot(x);
% hold on;
% plot(x_TBP,'m');
% title(['(b) Throw-away BP']);
% 
% subplot(4,1,3);
% plot(x);
% hold on;
% plot(x_DBP, 'g');
% title(['(c) Density BP']);
% 
% subplot(4,1,4);
% plot(x);
% hold on;
% plot(x_TDBP,'r');
% title(['(d) Throw-away Density BP']);
% 
% figure;
% subplot(4,1,1);
% plot(x);
% hold on;
% plot(x_OMP, 'r');
% title(['(a) OMP']);
% 
% subplot(4,1,2);
% plot(x);
% hold on;
% plot(x_DOMP,'m');
% title(['(b) Density OMP']);
% 
% subplot(4,1,3);
% plot(x);
% hold on;
% plot(x_TDOMP, 'g');
% title(['(c) Throw-away Density OMP']);
% 
% subplot(4,1,4);
% plot(x);
% hold on;
% plot(x_BCS,'r');
% title(['(d) Bayesian Compressive Sensing']);

disp(['BP: ||I_hat-I||/||I|| = ' num2str(E_BP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_BP) ', time = ' num2str(t_BP) ' secs']);
disp(['TBP: ||I_hat-I||/||I|| = ' num2str(E_TBP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_TBP) ', time = ' num2str(t_TBP) ' secs']);
disp(['DBP: ||I_hat-I||/||I|| = ' num2str(E_DBP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_DBP) ', time = ' num2str(t_DBP) ' secs']);
disp(['TDBP: ||I_hat-I||/||I|| = ' num2str(E_TDBP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_TDBP) ', time = ' num2str(t_TDBP) ' secs']);

disp(['OMP: ||I_hat-I||/||I|| = ' num2str(E_OMP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_OMP) ', time = ' num2str(t_OMP) ' secs']);
disp(['TOMP: ||I_hat-I||/||I|| = ' num2str(E_TOMP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_TOMP) ', time = ' num2str(t_TOMP) ' secs']);
disp(['DOMP: ||I_hat-I||/||I|| = ' num2str(E_DOMP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_DOMP) ', time = ' num2str(t_DOMP) ' secs']);
disp(['TDOMP: ||I_hat-I||/||I|| = ' num2str(E_TDOMP) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_TDOMP) ', time = ' num2str(t_TDOMP) ' secs']);
disp(['BCS: ||I_hat-I||/||I|| = ' num2str(E_BCS) ', ||y_hat - y|| / ||y|| = ' num2str(Ey_BCS) ', time = ' num2str(t_BCS) ' secs']);


