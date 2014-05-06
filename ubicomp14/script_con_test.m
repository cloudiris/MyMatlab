clear
clc

% s = load('con_09_11_seletcted.txt');
s = load('d:\document\codes\data\construction\c_test.txt');

g = 1:20;
lg = zeros(1, length(g));
pow = ones(1, length(g));

for i = length(g):-1:1
    lg(i) = length(unique(s(:,g(i))));
    pow(i) = prod(lg(i+1:end));
end

ind = s(:,g) * pow';