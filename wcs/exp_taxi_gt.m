clear
clc

workpath = 'd:\document\codes\matlab\wcs\';
seg = load('200seg.txt');

si1 = 401;
si2 = 428;

n = length(seg);
intvn = 26;

filename = ['200gt0517.txt';'200gt0518.txt';'200gt0519.txt';'200gt0521.txt';'200gt0522.txt';'200gt0523.txt';'200gt0524.txt'];
gt_all = zeros(intvn * 7, n);

for no = 1 : 7
    raw = load(['taxidata\' filename(no, :)]);
    gt_all((no - 1) * intvn + 1 : no * intvn, :) = raw;
end

f1 = fopen('200gt.txt', 'w');
for i = 1 : 7 * intvn
    for j = 1 : n
        fprintf(f1, '%f ', gt_all(i, j));
    end
    fprintf(f1, '\r\n');
end
fclose(f1);