% figure;
figure;
hold on;
xidx = 0.7:0.05:0.95;
interpm = 'spline';
if (exist('res_ca', 'var') && numel(res_ca) > 0)
    [~, in] = sort(res_ca(:,1));
    res_ca = res_ca(in, :);
    s_ca(1, :) = xidx;
    s_ca(2, :) = interp1(res_ca(:,6), res_ca(:,2), xidx, interpm);
    plot(xidx,s_ca(2,:), 'b-'); 
end

if (exist('res_pair', 'var') && numel(res_pair) > 0)
    [~, in] = sort(res_pair(:,1));
    res_pair = res_pair(in, :);
    s_pair(1, :) = xidx;
    s_pair(2, :) = interp1(res_pair(:,6), res_pair(:,2), xidx, interpm);
    plot(xidx,s_pair(2,:),'r--'); 
end

if (exist('res_tmc', 'var') && numel(res_tmc) > 0) 
    [~, in] = sort(res_tmc(:,1));
    res_tmc = res_tmc(in, :);
    s_tmc(1, :) = xidx;
    s_tmc(2, :) = interp1(res_tmc(:,6), res_tmc(:,2), xidx,interpm);
    plot(xidx,s_tmc(2,:),'k:.'); 
end

if (exist('res_random', 'var') && numel(res_random) > 0) 
    [~, in] = sort(res_random(:,1));
    res_random = res_random(in, :);
    s_random(1, :) = xidx;
    s_random(2, :) = interp1(res_random(:,6), res_random(:,2), xidx, interpm);
    plot(xidx,s_random(2,:),'m:'); 
end

if (exist('res_greedy', 'var') && numel(res_greedy) > 0) 
    [~, in] = sort(res_greedy(:,1));
    res_greedy = res_greedy(in, :);
    s_greedy(1, :) = xidx;
    s_greedy(2, :) = interp1(res_greedy(:,6), res_greedy(:,2), xidx, interpm);
    plot(xidx,s_greedy(2,:),'c--.'); 
end

if (exist('res_local', 'var') && numel(res_local) > 0) 
    [~, in] = sort(res_local(:,1));
    res_local = res_local(in, :);
    s_local(1, :) = xidx;
    s_local(2, :) = interp1(res_local(:,6), res_local(:,2), xidx, interpm);
    plot(xidx,s_local(2,:),'c--.'); 
end