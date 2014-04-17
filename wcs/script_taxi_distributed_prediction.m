clear
clc

NDEBUG = 1;

workpath = 'd:\document\codes\matlab\wcs\';
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
combine = 1;
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
mfrac_min = 0.5;
mfrac_gap = 0.05;
mfrac_max = 0.5;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracs = [9.70E-01	9.37E-01	9.00E-01	8.59E-01	8.13E-01	7.61E-01	7.04E-01	6.39E-01	5.67E-01	4.86E-01	3.98E-01	3.04E-01	2.09E-01	1.20E-01	5.10E-02	1.23E-02	9.31E-04	3.91E-06	6.99E-13];
% mfracs = [3.91E-02,5.22E-03];
mfracn = length(mfracs);

times = 1;
single_runs = mfracn * times;

samplings = 1;
alpha = 2.5;
samplingn = length(samplings);

runs = samplingn * mfracn;

exp_taxi_distributed_prediction;

save('res_taxi_prediction.txt', 'res', '-ascii');