clear
clc

s = load('d:\document\codes\data\pedestrian\ped.txt');
sw = [43.5917, -79.6325];

for i = 1 : length(s)
    s(i,1:2) = [ceil(distance(sw(1), sw(2), sw(1), s(i,2), 'degree')/360*40000), ...
                ceil(distance(sw(1), sw(2), s(i,1), sw(2), 'degree')/360*40000)];
    
end
