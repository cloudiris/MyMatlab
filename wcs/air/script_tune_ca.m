%% tune CA
% clear
clc

% pm25i = load('pm25i.txt');
dim = size(pm25i);
n = dim(2);
d = sqrt(n);
k = 50;
b = base_dct4(n);
snapn = 1;
%snaps = sort(randsample(dim(1), snapn));
snapstart = 3666;
snaps = snapstart : snapstart + snapn - 1;
fprintf('data loaded: %dx%d\n', dim(1), dim(2));

cost = load('D:\Document\Codes\Matlab\wcs\3gcost.32.32.txt');
cost = cost(1:d, 1:d);
costv = hilbertcurve(cost);
sumcost = sum(costv);
fprintf('cost loaded\n');

alp_min = 0;
alp_gap = 0.2;
alp_max = 5;
alps = alp_min:alp_gap:alp_max;
alpn = length(alps);

m_min = 50;
m_max = n;
m_gap = round(n/10);
ms = m_min : m_gap : m_max;
mn = length(ms);

denoise = 0;
tolerance = 25;
times = 10;
totaltimes = times * alpn * snapn * mn;

res = zeros(totaltimes, 5); % m|alp|cost|acc|snap

loopcounter = 0;

%% main loop

for alpi = 1 : alpn;
	alp = alps(alpi);
	prob = costv.^(-alp) / sum(costv.^(-alp));
	
	for mi = 1 : mn
		m = ms(mi);
		
		for timei = 1 : times
			%%sampling
			samplex = randsamplewtr(n, m, prob);
			costratio = sum(costv(samplex)) / sumcost;
			
			%%recovery
			for snapi = 1 : snapn
				loopcounter = loopcounter + 1;
				snap = snaps(snapi);
				s = pm25i(snap, :)';
				re = my_csf(s(samplex), samplex, b, k, 'OMP', denoise);
				re(re < 0) = 0;
                acc = sum(abs(re-s) <= tolerance) / n;
				res(loopcounter, :) = [m, alp, costratio, acc, snap];
			end
		end
		fprintf('m(%d)|alp(%f) is over\n', m, alp);
	end
	fprintf('alp(%f) is over\n', alp);
end

imethod = 'linear';
coststops = 0:0.01:1;
sres_acc = [0 coststops; zeros(alpn, 1+length(coststops))];
accstops = 0:0.01:1;
sres_cost = [0 accstops; zeros(alpn, 1+length(accstops))];
for alpi = 1:alpn
    alp = alps(alpi);
    mres = zeros(mn, 2);
    for mi = 1:mn
        m = ms(mi);
        mres(mi, :) = mean(res(res(:,1)==m & res(:,2)==alp, [3 4]));
    end
    sres_acc(alpi+1, 1) = alp;
    sres_acc(alpi+1, 2 : 1+ length(coststops)) = interp1([0 mres(:, 1)' 1], [0 mres(:,2)' 1], coststops, imethod);
    
    [sacc, idx] = sort(mres(:,2));
    [usacc, uidx] = unique(sacc);
    idx = idx(uidx);
    sres_cost(alpi+1, 1) = alp;
    sres_cost(alpi+1, 2 : 1+ length(accstops)) = interp1([0 mres(idx, 2)' 1], [0 mres(idx,1)' 1], accstops, imethod);
end

figure
