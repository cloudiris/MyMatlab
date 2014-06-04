function bases = train_pm25(y, train_per, ks)

dim = size(y);

n = dim(1); %measurment dimension = 1024
N = dim(2); %# samples = 5000

trainstop = floor(N * train_per);

data = y(:, 1:trainstop);

k_size = size(ks, 1);
bases = zeros(k_size, n, n);

for i = 1:k_size
    tb = orth(simple_ksvd(data, n, ks(i), 1));
    b(i, 1:size(tb,1),1:size(tb,2)) = tb;
    fprintf('%d is over, rank=%d\n', k, rank(tb));
end