clear
clc

s = load('d:\document\codes\data\nonemployment\non.txt');
loc = unique(s(:,1)*1000 + s(:,2));
type = unique(s(:,3));

ld1 = length(type);
ld2 = length(loc);

summary1 = zeros(ld1, ld2);
summary2 = summary1;
count = summary1;

for i = 1 : length(s)
    i1 = find(type==s(i,3));
    i2 = find(loc == s(i,1)*1000+s(i,2));
    summary1(i1,i2) = summary1(i1,i2) + s(i,4);
    summary2(i1,i2) = summary2(i1,i2) + s(i,5);
    count(i1,i2) = count(i1,i2) + 1;
    fprintf('%d\n', i);
end

idx = find(count(:)>0);
summary1(idx) = summary1(idx)./count(idx);
summary2(idx) = summary2(idx)./count(idx);

rsum = sum(count, 2);
csum = sum(count, 1);
[~,idx1] = sort(rsum, 'descend');
[~,idx2] = sort(csum, 'descend');
idx1 = idx1(randperm(150));
idx2 = idx2(randperm(1000));

s1 = summary1(idx1, idx2);
s2 = summary2(idx1, idx2);

train1 = s1(:,1:700); test1 = s1(:,701:1000);
train2 = s2(:,1:700); test2 = s2(:,701:1000);
