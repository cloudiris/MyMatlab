clear
clc

NDEBUG = 0;
%% pretreatment
N = 10000;
low = 0;
high = 100;

% cost = low + (high - low) .* rand(1, N);
% cost = normrnd(3,1, 1,N);
cost = exprnd(1, 1, N);

cost = cost(find(cost > low));
cost = cost(find(cost < high));
N = numel(cost);
sumcost = sum(cost);

%% parameter
recovery = 'CS';
mfrac_min = 0.1;
mfrac_gap = 0.05;
mfrac_max = 0.9;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracn = length(mfracs);

times = 50;
single_runs = mfracn * times;

samplings = [0,1,1.1,1.2,1.3,1.4,2,2.3,2.4];
samplingn = length(samplings);

runs = samplingn * single_runs;

exp_analysis;

res_random = res(find(res(:,1)==0), 2:size(res,2));
res_ca1 = res(find(res(:,1)==1), 2:size(res,2));
res_cadot1 = res(find(res(:,1)==1.1), 2:size(res,2));
res_cadot5 = res(find(res(:,1)==1.2), 2:size(res,2));
res_ca2 = res(find(res(:,1)==1.3), 2:size(res,2));
res_ca5 = res(find(res(:,1)==1.4), 2:size(res,2));
res_2group = res(find(res(:,1)==2), 2:size(res,2));
res_3group = res(find(res(:,1)==2.3), 2:size(res,2));
res_4group = res(find(res(:,1)==2.4), 2:size(res,2));

d2 = size(res_random, 2);

f = fopen('res_analysis.txt', 'w');
sf = fopen('res_analysis_s.txt', 'w');
for samplingi = 1 : samplingn
    sampling = samplings(samplingi);
    if (sampling == 0) 
        output_res = res_random;
        fprintf(f, 'res_random\n');
        s_r = intervaled(res_random(:, 2:d2), mfracn);
        s_output = s_r;
        fprintf(sf, 'random\n');
    elseif (sampling == 1)
        output_res = res_ca1;
        fprintf(f, 'res_ca alpha=1\n');
        s_ca1 = intervaled(res_ca1(:, 2:d2), mfracn);
        s_output = s_ca1;
        fprintf(sf, 'ca 1\n');
    elseif (sampling == 1.1)
        output_res = res_cadot1;
        fprintf(f, 'res_ca alpha=0.1\n');
        s_cadot1 = intervaled(res_cadot1(:, 2:d2), mfracn);
        s_output = s_cadot1;
        fprintf(sf, 'ca 0.1\n');
    elseif (sampling == 1.2)
        output_res = res_cadot5;
        fprintf(f, 'res_ca alpha=0.5\n');
        s_cadot5 = intervaled(res_cadot5(:, 2:d2), mfracn);
        s_output = s_cadot5;
        fprintf(sf, 'ca 0.5\n');
    elseif (sampling == 1.3)
        output_res = res_ca2;
        fprintf(f, 'res_ca alpha=2\n');
        s_ca2 = intervaled(res_ca2(:, 2:d2), mfracn);
        s_output = s_ca2;
        fprintf(sf, 'ca 2\n');
    elseif (sampling == 1.4)
        output_res = res_ca5;
        fprintf(f, 'res_ca alpha=5\n');
        s_ca5 = intervaled(res_ca5(:, 2:d2), mfracn);
        s_output = s_ca5;
        fprintf(sf, 'ca 5\n');
    elseif (sampling == 2)
        output_res = res_2group;
        fprintf(f, 'res_2group\n');
        s_2group = intervaled(res_2group(:, 2:d2), mfracn);
        s_output = s_2group;
        fprintf(sf, '2group\n');
    elseif (sampling == 2.3)
        output_res = res_3group;
        fprintf(f, 'res_3group\n');
        s_3group = intervaled(res_3group(:, 2:d2), mfracn);
        s_output = s_3group;
        fprintf(sf, '3group\n');
    elseif (sampling == 2.4)
        output_res = res_4group;
        fprintf(f, 'res_4group\n');
        s_4group = intervaled(res_4group(:, 2:d2), mfracn);
        s_output = s_4group;
        fprintf(sf, '4group\n');
    else
        continue;
    end
    [rc,cc] = size(output_res);
    for i = 1 : rc
        for j = 1 : cc
            fprintf(f, '%f ', output_res(i,j));
        end
        fprintf(f, '\n');
    end
    
    [rc,cc] = size(s_output);
    for i = 1 : rc
        for j = 1 : cc
            fprintf(sf, '%f ', s_output(i,j));
        end
        fprintf(sf, '\n');
    end
end
fclose(f);
fclose(sf);