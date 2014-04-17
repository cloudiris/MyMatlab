clear
clc

datapath = 'd:\document\codes\data\';
workpath = 'd:\document\codes\matlab\wcs\';

position = load([datapath 'newposition.txt']);
gt = load([datapath 'gt20090525.txt']);
plot(position(:,1), position(:, 2), 'o');

