%% s

samplings = [0 1 4 5];
samplingn = length(samplings);

accstops = 0:0.01:1;
sres_target = zeros(1,length(accstops));

imethod = 'linear';
for samplingi = 1 : samplingn
    sampling = samplings(samplingi);
    
    if(sampling == 0) tres = res_random;
    elseif (sampling == 1) tres = res_ca;
    elseif (sampling == 4) tres = res_tmc;
    elseif (sampling == 5) tres = res_greedy; end
    [sacc, idx] = sort(tres(:,2));
    [usacc, uidx] = unique(sacc);
    idx = idx(uidx);
    sres_target = interp1([0 tres(idx, 2)' 1], [0 tres(idx, 6)' 1], accstops, imethod);
    
    if(sampling == 0) sres_random = sres_target;
    elseif (sampling == 1) sres_ca = sres_target;
    elseif (sampling == 4) sres_tmc = sres_target;
    elseif (sampling == 5) sres_greedy = sres_target; end
end

sres_all = [sres_random; sres_ca; sres_tmc; sres_greedy];
idx = [81 86 91 96];
figure
bar(sres_all(:, idx)');