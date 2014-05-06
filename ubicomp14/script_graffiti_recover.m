clear
clc

load('graffiti_base_month.mat');
test = load('graffiti_test_month.txt');

% k = 6;
n = 100;
m_min = 10;
m_gap = 10;
m_max = 90;
ms = m_min : m_gap : m_max;
mn = length(ms);

methods = [1 2 3 4]; % 1 cs; 2 fourier; 3 linear; 4 spline
methodn = length(methods);

times = 10;

dctb = base_dct4(n);
raw = zeros(n, size(test,2), mn*methodn*times);
config = zeros(mn * times, 2); % m|method
tres = zeros(size(test,2), 3);
res = zeros(mn, 4);
innerc = 0;

for timei = 1 : times
    
    for mi = 1 : mn
        m = ms(mi); 
        k = round(m/log2(n));
    
        samplex = randperm(n, m);
        ssamplex = sort(samplex);
    
        for methodi = 1 : methodn
            method = methods(methodi);
            innerc = innerc + 1;
            config(innerc, :) = [m, method];
            
            for si = 1 : size(test,2)
                s = test(:,si);
%                 tolv = s .* tol;
                
                if (method ==1)
                    re = my_csf(s(samplex),samplex,b(:,:,k),k,'OMP',1);
                elseif (method == 2)
                    re = my_csf(s(samplex),samplex,dctb,k,'OMP',1);
                elseif (method == 3)
                    re = interp1(ssamplex, s(ssamplex), 1:n, 'linear', 'extrap')';
                elseif (method == 4)
                    re = interp1(ssamplex, s(ssamplex), 1:n, 'spline', 'extrap')';
                end

                re(re<0) = 0;
%                 re = round(re);
                raw(:,si,innerc) = re;
%                 acc1 = snr(s, re1(:,innerc));
%                 acc2 = snr(s, re2(:,innerc));
%                 acc3 = snr(s, re3(:,innerc));
% 
%                 acc1 = sum(abs(s-re1(:,innerc))<=tolv)/n;
%                 acc2 = sum(abs(s-re2(:,innerc))<=tolv)/n;
%                 acc3 = sum(abs(s-re3(:,innerc))<=tolv)/n;
% 
%                 acc1 = 1-sum(abs(s-re1(:,innerc)))/sum(s);
%                 acc2 = 1-sum(abs(s-re2(:,innerc)))/sum(s);
%                 acc3 = 1-sum(abs(s-re3(:,innerc)))/sum(s);
% 
%                 acc1 = 1-norm(abs(s-re1(:,innerc)))/norm(s);
%                 acc2 = 1-norm(abs(s-re2(:,innerc)))/norm(s);
%                 acc3 = 1-norm(abs(s-re3(:,innerc)))/norm(s);
%                 tres(innerc,1:3) = [acc1,acc2,acc3];
            end
        end
    end
%     res(mi, :) = [m, mean(tres)];
end

re_by_method = zeros(n, size(test,2), mn, methodn);
for mi = 1 : mn
    m = ms(mi);
    for methodi = 1 : methodn
        method = methods(methodi);
        idx = find(config(:,1)==m & config(:,2)==method);
        re_by_method(:,:,mi,methodi) = mean(raw(:,:,idx), 3);
    end
end

%% separate
tol = 0.4;
sep_res = zeros(mn, methodn);
for mi = 1 : mn
    for methodi = 1 : methodn
        for si = 1 : size(test,2)
            s = test(:,si);
            tolv = s.*tol;
%             sep_res(mi, methodi) = sep_res(mi, methodi) + snr(s, re_by_method(:,si,mi,methodi));
            sep_res(mi, methodi) = sep_res(mi, methodi) + sum(abs(s-re_by_method(:,si,mi,methodi))<tolv)/n;
        end
        sep_res(mi, methodi) = sep_res(mi, methodi) / size(test,2);
    end
end
figure
plot(ms, sep_res);
legend('CS w/ trained base', 'CS w/ Fourier base', 'Linear interp.', 'Spline interp.');
%% total number
% totaltruth = mean(test, 2);
% tolv = totaltruth.*tol;
% total = zeros(n, mn, methodn);
% total_res = zeros(mn, methodn);
% for mi = 1 : mn
%     m = ms(mi);
%     for methodi = 1 : methodn
%         method = methods(methodi);
%         total(:,mi,methodi) = mean(re_by_method(:,:,mi,methodi), 2);
% %         total_res(mi, methodi) = snr(totaltruth, total(:,mi, methodi));
%         total_res(mi, methodi) = sum(abs(totaltruth-total(:,mi, methodi))<tolv);
%     end
% end
% figure;
% plot(ms, res(:,2:4));
% legend('cs','fourier','interp');