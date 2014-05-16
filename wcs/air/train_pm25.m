function bases = train_pm25(y, ks)

fprintf('training start\n');
dim = size(y);

n = dim(1); %measurment dimension = 100
N = dim(2); %# samples = 400
K = floor(n * 1.5);

k_size = size(ks, 2);
bases = zeros(k_size, n, K);

for i = 1:k_size
    k = ks(i);
    tb = orth(simple_ksvd(y, K, k, 1));
    tb(1:10, 1:10)
    bases(1:size(tb,1),1:size(tb,2),i) = tb;
    bases(1:10, 1:10, i)
    fprintf('%d is over, rank=%d\n', k, rank(tb));
end