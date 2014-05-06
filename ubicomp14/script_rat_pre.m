clear
clc

s = load('d:\document\codes\data\NYC-rat\east_rat.txt');
summary = zeros(40);

for i = 1 : length(s)
    x = ceil((s(i,1)-970000)/2500);
    y = ceil((s(i,2)-140000)/3500);
    summary(x,y) = summary(x,y) + 1;
    fprintf('%d\n',i);
end