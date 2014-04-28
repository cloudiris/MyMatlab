clear
clc

NDEBUG = 0;

if (ispc()) workpath = 'd:\document\codes\matlab\wcs\';
else workpath = '~/Documents/Codes/Matlab/wcs/'; end

week = 'w2';
seg = load('200seg.txt');
n = length(seg);
gt = load(['200gti_' week '.txt']);
[intvn, ~] = size(gt);
% intvn = 26;
% positioncost = load('200positioncost_conn.txt');
% positioncost = load('gpscost.txt');
% positioncost = load('3gcost.txt');
positioncost = load('gps+3g.txt');
% positioncost = ones(101,107);
dim = size(positioncost);
% timecost = load(['battercost_li_' week '.txt']);
timecost = load(['battercost_ex_' week '.txt']);
% timecost = load(['battercost_lo_' week '.txt']);
% timecost = load(['30Mtimecost_' week '.txt']);
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
mfrac_min = 0.05;
mfrac_gap = 0.1;
mfrac_max = 1;
mfracs = mfrac_min:mfrac_gap:mfrac_max;
mfracn = length(mfracs);

times = 1;
single_runs = mfracn * times;

samplings = [0 1 2 3 4 5]; %0 random, 1 ca, 2 pair, 3 pw-d, 4 tmc, 5 greedy
alp = 2.7;
dis_thr = 1; % in km
% k = 2;
samplingn = length(samplings);

runs = samplingn * single_runs;

exp_taxi_ca;

% afrac|acc|err|coverage|totalcost|ratiocost|time
res_random = res(find(res(:,1)==0), 2:8);
res_ca = res(find(res(:,1)==1), 2:8);
res_pair = res(find(res(:,1)==2), 2:8);
res_pw_d = res(find(res(:,1)==2.1), 2:8);
res_tmc = res(find(res(:,1)==4), 2:8);
res_greedy = res(find(res(:,1)==5), 2:8);
if (size(res_random,1) == mfracn && size(res_random, 1) == size(res_greedy, 1)) 
    res_optimal = res_random;
    res_optimal(:,[5 6]) = res_greedy(:, [5 6]);
end

imethod = 'linear';
accstops = 0:0.01:1;
sres_target = zeros(length(accstops), 1);

for samplingi = 1 : samplingn
    sampling = samplings(samplingi);
    if (sampling == 0) tres = res_random;
    elseif (sampling == 1) tres = res_ca;
    elseif (sampling == 2) tres = res_pair;
    elseif (sampling == 2.1) tres = res_pw_d;
    elseif (sampling == 4) tres = res_tmc;
    elseif (sampling == 5) tres = res_greedy;
    end
    [sacc, idx] = sort(tres(:,3));
    [usacc, uidx] = unique(sacc);
    idx = idx(uidx);
    sres_target = interp1([0 tres(idx,2)' 1], [0 tres(idx,6)' 1], accstops, imethod)';
    if (sampling == 0) sres_random = sres_target;
    elseif (sampling == 1) sres_ca = sres_target;
    elseif (sampling == 2) sres_pair = sres_target;
    elseif (sampling == 2.1) sres_pw_d = sres_target;
    elseif (sampling == 4) sres_tmc = sres_target;
    elseif (sampling == 5) sres_greedy = sres_target;
    end
end

% interpolate optimal
if (exist('res_optimal','var'))
    tres = res_optimal;
    [sacc, idx] = sort(tres(:,3));
    [usacc, uidx] = unique(sacc);
    idx = idx(uidx);
    sres_target = interp1([0 tres(idx,2)' 1], [0 tres(idx,6)' 1], accstops, imethod)';
    sres_optimal = sres_target;
end

figure
hold on
if (exist('sres_random','var')) plot(accstops,sres_random); end
if (exist('sres_pair','var')) plot(accstops,sres_pair, 'r'); end
if (exist('sres_pw_d','var')) plot(accstops,sres_pw_d, 'c'); end
if (exist('sres_ca','var')) plot(accstops,sres_ca, 'm'); end
if (exist('sres_greedy','var')) plot(accstops,sres_greedy, 'k'); end
if (exist('sres_llc','var')) plot(accstops,sres_llc, 'g'); end
if (exist('sres_optimal','var')) plot(accstops, sres_optimal, 'y'); end

sres_all = [];
if (exist('sres_random','var')) sres_all = [sres_all sres_random]; end
if (exist('sres_pair','var')) sres_all = [sres_all sres_pair]; end
if (exist('sres_pw_d','var')) sres_all = [sres_all sres_pw_d]; end
if (exist('sres_ca','var')) sres_all = [sres_all sres_ca]; end
if (exist('sres_tmc','var')) sres_all = [sres_all sres_tmc]; end
if (exist('sres_greedy','var')) sres_all = [sres_all sres_greedy]; end
if (exist('sres_llc','var')) sres_all = [sres_all sres_llc]; end
if (exist('sres_optimal','var')) sres_all = [sres_all sres_optimal]; end
idx = [76 81 86 91 96];
figure;
bar(sres_all(idx, :));

% f = fopen('res_taxi.txt', 'w');
% for samplingi = 1 : samplingn
%     sampling = samplings(samplingi);
%     if (sampling == 0) 
%         output_res = res_random;
%         fprintf(f, 'res_random\n');
%     elseif (sampling == 1)
%         output_res = res_ca;
%         fprintf(f, 'res_ca\n');
%     elseif (sampling == 2)
%         output_res = res_pair;
%         fprintf(f, 'res_pair\n');
%     elseif (sampling == 4)
%         output_res = res_tmc;
%         fprintf(f, 'res_tmc\n');
%     elseif (sampling == 5)
%         output_res = res_greedy;
%         fprintf(f, 'res_greedy\n');
%     else
%         continue;
%     end
%     [rc,~] = size(output_res);
%     for i = 1 : rc
%         for j = 1 : 7
%             fprintf(f, '%f ', output_res(i,j));
%         end
%         fprintf(f, '\n');
%     end
% end
% fclose(f);