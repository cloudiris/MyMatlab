clear
clc

set = 2;

train = load(['non_train' int2str(set) '.txt']);
test = load(['non_test' int2str(set) '.txt']);

N = 150;
k = 20;
b = zeros(N, N, k);
for i = 1 : k
    tb = orth(simple_ksvd(train, N, i, 1));
    b(1:size(tb,1), 1:size(tb,2), i) = tb;
    fprintf('%d over\n', i);
end

save(['non_base' int2str(set) '.mat'],'b','-mat');