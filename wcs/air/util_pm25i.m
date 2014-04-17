clear
clc

inc = [2 3 4 6 7 8 9 10 11 12 13 14 15 16 18 20 21 22];
loc = load('station_coor.txt');
s = load('pm25.txt');
s = s(:, inc);
for i = 1 : size(s,2)
    idx = find(s(:,i)~=-1);
    s(:,i) = interp1(idx, s(idx, i), 1:size(s,1), 'linear', 'extrap');
end

for i = 1 : size(s, 1)
    s(i, 19) = s(i, randi(18));
end

pm25i = zeros(size(s, 1), 1024);
[X,Y] = meshgrid(1:32, 1:32);
for i = 1 : size(s, 1)
    pm25i(i, :) = hilbertcurve(griddata(loc(:,1), loc(:,2), s(i, :), X, Y, 'v4'))';
    fprintf('%d\n', i);
end

pm25i(find(pm25i<0)) = -1*pm25i(find(pm25i<0));