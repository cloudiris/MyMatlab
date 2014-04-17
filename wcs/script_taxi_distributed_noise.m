clear
clc

NDEBUG = 0;

workpath = 'e:\document\codes\matlab\wcs\';
week = 'w2';
seg = load('200seg.txt');
n = length(seg);
gt = load(['200gti_' week '.txt']);
[intvn, ~] = size(gt);
% intvn = 26;
% positioncost = load('200positioncost_conn.txt');
positioncost = load('gpscost.txt');
% positioncost = ones(101,107);
dim = size(positioncost);
% timecost = load(['battercost_li_' week '.txt']);
timecost = load(['battercost_ex_' week '.txt']);
% timecost = load(['battercost_lo_' week '.txt']);
% timecost = ones(9000, intvn);
combine = 1; % 0 product, 1 sum;
recnum = load(['30Mrecnum_' week '.txt']);
tN = sum(recnum(1:intvn));
cabcount = load(['cabcount_' week '.txt']);
cabID = find(cabcount  > 0);
cabnum = numel(cabID);

segmapping = zeros(dim(1), dim(2));
for pt = 1 : n
    segmapping(seg(pt,1) + 1, seg(pt,2) + 1) = pt;
end

recovery = 'CS';
intvbatch = intvn;

mfracs = [9.00E-01	8.59E-01	8.13E-01	7.61E-01	7.04E-01	6.39E-01	5.67E-01	4.86E-01];
% mfracs = 0.3:0.05:0.50;
mfracn = length(mfracs);

times = 1;
single_runs = mfracn * times;

noises = [-15 -10 -5 -1];
noisen = numel(noises);

alp = 1;
kk = 1;
sampling = 1;

runs = noisen * single_runs;

exp_taxi_distributed_noise;

save('res_taxi_distributed_noise.txt','res','-ascii');