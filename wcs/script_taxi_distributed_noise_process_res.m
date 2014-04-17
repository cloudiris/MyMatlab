interpm = 'linear';

acctointerp = 0.8:0.05:0.9;
noisetointerp = 0.1:0.1:1.0;

res_noise = zeros(numel(noisetointerp), numel(acctointerp));

tnoise = unique(res(:,1));
for i = 1 : numel(tnoise)
    no = tnoise(i);
    res_tmp = res(find(res(:,1)==tnoise(i)), 2:8);
    clear icost;
    icost = interp1(res_tmp(:,2)', res_tmp(:,6)', acctointerp, interpm, 'extrap');
    res_noise(i, 1:numel(acctointerp)) = icost;
end
% figure;
% plot(res_noise())
% 
% for i = 1 : numel(acctointerp)
%     clear icost;
%     icost = interp1(tnoise',res_noise(1:numel(tnoise), i), noisetointerp, interpm, 'extrap');
%     res_noise(:,i) = icost;
% end
% figure;
% plot(noisetointerp, res_noise());
% 
% res_noise = [noisetointerp' res_noise];
save('res_taxi_noise.txt', 'res_noise', '-ascii');