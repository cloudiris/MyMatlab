function [y_recover x_recover] = my_wbp_easy(y, sigma, A, eps)

A_full = A;

mask = find(~isnan(y));
y = y(mask);
sigma = sigma(mask);
A = A(mask, :);

Ab = inv(diag(sigma .^ 2)) * A;
yb = inv(diag(sigma .^ 2)) * y;
epsilon = sqrt(size(y(mask), 1)) * eps;
x0 = Ab'*inv(Ab*Ab')*yb;
x_DBP = l1qc_logbarrier(x0, Ab, [], yb, epsilon, 1e-3);
y_DBP = A_full * x_DBP; 

fprintf(1,'DBP number of nonzero weights: %d %f\n',sum(x_DBP~=0), sigma(1));
x_recover = x_DBP;
y_recover = y_DBP;