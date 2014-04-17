% clear
% clc
n = 10000;
m = 1000;

w = abs(markovc(n))'*10;
% w = sort(w);
% w = w / sum(w);

iter = 1000;
counter = zeros(1, n);
add = ones(1, m);
tic;
for i = 1 : iter
%     idx = randsamplewtr(n, m, w);
    key = (rand(1, n)).^(1./w);
    [~, big] = sort(key, 'descend');
    idx = big(1:m);
    counter(idx) = counter(idx) + add;
end
toc;
counter = counter / iter / m;

figure;
plot(counter);
hold on
plot(w/sum(w), 'r');