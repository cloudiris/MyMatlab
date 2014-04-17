clear
% clc
NDEBUG = 0;
%% pretreatment
workpath = 'd:\document\codes\matlab\wcs\';
datapath = 'd:\document\codes\data\';
data = 'temp.grid.64.64.txt';
s2 = load([datapath data]);
s2(:,1:32) = -58;
s = hilbertcurve(s2)';
N = numel(s);
dim = size(s2);

% positioncost = load('exp1.txt');
% positioncost = load('exp2.txt');
% positioncost = load('exp3.txt');
positioncost = load('local_example.txt');
% positioncost = load('local_example3.txt');

cost = hilbertcurve(positioncost);
sumcost = sum(cost);
b = base_dct4(N);

recovery = 'CS';

mfrac_min = 0.05;
mfrac_gap = 0.05;
mfrac_max = 0.70;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracn = length(mfracs);

times = 5;
single_runs = mfracn * times;

samplings = [1];
samplingn = length(samplings);

runs = samplingn * mfracn;

exp_2d;

fprintf('afrac=%f, acc=%f, cost=%f\n', res(1,[2 5 7]));

% res_random = res(find(res(:,1)==0), 2:8);
res_ca = res(find(res(:,1)==1), 2:8);
% res_pair = res(find(res(:,1)==2), 2:8);
res_greedy = res(find(res(:,1)==5), 2:8);
res_local = res(find(res(:,1)==6), 2:8);

f = fopen('res_2d.txt', 'w');
for samplingi = 1 : samplingn
    sampling = samplings(samplingi);
    if (sampling == 0) 
        output_res = res_random;
        fprintf(f, 'res_random\n');
    elseif (sampling == 1)
        output_res = res_ca;
        fprintf(f, 'res_ca\n');
    elseif (sampling == 2)
        output_res = res_pair;
        fprintf(f, 'res_pair\n');
    elseif (sampling == 5)
        output_res = res_greedy;
        fprintf(f, 'res_greedy\n');
    else
        continue;
    end
    [rc,~] = size(output_res);
    for i = 1 : rc
        for j = 1 : 7
            fprintf(f, '%f ', output_res(i,j));
        end
        fprintf(f, '\n');
    end
end
fclose(f);
figure;
hold on;
plot(res_ca(:,6),res_ca(:,4), 'b.');
plot(res_local(:,6),res_local(:,4), 'r.');
% plot(res_random(:,6),res_random(:,4), 'm--');
plot(res_greedy(:,6),res_greedy(:,4), 'k--');