clear
clc

% s = load('d:\document\codes\data\construction\c_all.txt');
% 
% figure
% for i = 1 : 30
%     subplot(5,6,i);
%     hist(s(:,i));
% end
% 
% figure
% for i = 1 : 30
%     subplot(5,6,i);
%     hist(s(:,30 + i));
% end

s = load('con_09_11_seletcted.txt');

g1 = [3 5 6 7 8 9];
g2 = [1 2 4 10];

lg1 = zeros(1, length(g1)); % length of dimensions in g1
pow1 = ones(1, length(g1)); % power base of dimensions in g1

lg2 = zeros(1, length(g2));
pow2 = ones(1, length(g2));

for i = length(g1) : -1 : 1
    lg1(i) = length(unique(s(:,g1(i))));
    pow1(i) = prod(lg1(i+1:end));
end
for i = length(g2) : -1 : 1
    lg2(i) = length(unique(s(:,g2(i))));
    pow2(i) = prod(lg2(i+1:end));
end

summary = zeros(prod(lg1), prod(lg2));
count = summary;

for i = 1 : length(s)
    c1 = 1 + (s(i,g1)-1)*pow1';
    c2 = 1 + (s(i,g2)-1)*pow2';
    summary(c1, c2) = summary(c1, c2) + s(i,12);
    count(c1, c2) = count(c1, c2) + 1;
end

idx = find(count(:)>0);
summary(idx) = summary(idx)./count(idx);

rsum = sum(count, 2);
csum = sum(count, 1);

[~, idx1] = sort(rsum, 'descend');
[~, idx2] = sort(csum, 'descend');

ss = summary(idx1(1:20), idx2(1:20));