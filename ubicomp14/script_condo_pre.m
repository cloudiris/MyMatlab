clear
clc 

year = 2008;
s = load(['d:\document\codes\data\nyc-condo\' int2str(year) '.txt']);

idx = find(s(:,3)~=0);
s = s(idx, :);
s(s(:,2) == 0, 2) = 10;

for i = 1 : length(s)
    s(i,3) = floor((year - s(i,3))/10);
    if (s(i,3) > 9) s(i,3) = 9; end
    
    s(i, 4) = floor((s(i,4) - 300)/300);
    if (s(i,4) < 0) s(i,4) = 0; end
    if (s(i,4) > 9) s(i,4) = 9; end
end

[~, idx] = sort(s(:,1));
s = s(idx, :);

d0 = unique(s(:,1));
d1 = unique(s(:,2));
ld0 = length(d0);
ld1 = length(d1);
ld2 = 10;
ld3 = 10;

summary = zeros(ld1*ld2*ld3, ld0);
count = summary;
for i = 1 : length(s)
    c0 = find(d0==s(i,1));
    c1 = find(d1==s(i,2)) - 1;
    c2 = s(i,3);
    c3 = s(i,4);
    cc = c1*(ld2*ld3) + c2*ld3 + c3 + 1;
    summary(cc, c0) = summary(cc, c0) + s(i,5);
    count(cc, c0) = count(cc, c0) + 1;
    fprintf('%d,%d\n', cc, c0);
end

idx = find(count(:)>0);
summary(idx) = summary(idx)./count(idx);

% save(['s' int2str(year) '.txt'], 's', '-ascii');