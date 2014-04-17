% clear
% clc

slpath = 'd:\document\codes\data\';
data = 'temp.grid.64.64';
datafile = [slpath data '.txt'];
s2 = load(datafile);
s = hilbertcurve(s2)';

si = size(s2);
n = si(1) * si(2);
k = floor(0.01 * n);
base = base_dct4(n);

m = 200;

recovery_method = 'OMP';

cost_distribution = 'station-uniform';
cost = generatecost(si(1), si(2), cost_distribution);
costv = hilbertcurve(cost);

%% finding optimal
%%% by annealing
% iteration = 1000;
% proc_annealing_2d;

%%% by file
optfile = [slpath 'optimal.data(' data ').cost(' cost_distribution ').mat'];
globaloptimal = load(optfile);

%% sampling parameter setup

times = 20;
sampling_method = 'ca';

samplex = zeros(1, m);
samplep = zeros(m, 2);

%% pick samples

if (sampling_method == 3 || sampling_method == 4) %cost-aware (plus random)
    pa = 0;
    proc_ca;
elseif (sampling_method == 5) % 2-D local least cost
    proc_2dllc;
    
elseif (sampling_method == 6) % 1-D LLC
    proc_1dllc;
    
elseif (sampling_method == 7) % evo w/ Info
    proc_evoinfo;
end

% s = hilbertcurve(rotm(s2', -90));
re = my_csf(s(samplex), samplex, base, k, recovery_method, 1);
% re2 = rotm(antihc(re), 90)';
totalcost = sum(costv(samplex));
costratio = totalcost / allcost;
cacc = acc100(s,re);
cr = cacc / costratio;

fprintf('Accracy = %.3f%%, Cost = %.3f%%, APC = %f\nOptimal = %f, Ratio = %f\n', cacc * 100, costratio * 100, cr, globaloptimal, cr / globaloptimal);