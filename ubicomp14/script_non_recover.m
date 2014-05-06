clear
clc

set = 1;
load(['non_base' num2str(set) '.mat']);
test = load(['non_test' num2str(set) '.txt']);

% k = 6;
n = 150;
m_min = round(0.1*n);
m_gap = round(0.1*n);
m_max = round(0.9*n);
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

            end
        end
    end
%     res(mi, :) = [m, mean(tres)];
end

tol = 0.3;

%% averaged acc by method
acc_by_method = zeros(mn, methodn);
for mi = 1 : mn
    m = ms(mi);
    for methodi = 1 : methodn
        method = methods(methodi);
        idx = find(config(:,1)==m & config(:,2)==method);
        acc = 0;
        for ini = 1 : length(idx)
            for si = 1 : size(test,2)
                s = test(:,si);
                tolv = s.*tol;
    %             sep_res(mi, methodi) = sep_res(mi, methodi) + snr(s, re_by_method(:,si,mi,methodi));
                acc = acc + sum(abs(s-raw(:,si,idx(ini)))<tolv)/n;
%                 acc = acc + (1-norm(s-raw(:,si,idx(ini)),1)/norm(s,1));
            end
        end
        acc = acc/length(idx)/size(test,2);
        acc_by_method(mi, methodi) = acc;
    end
end
figure;
plot(ms/n, acc_by_method);
legend('CCS', 'CS', 'Linear interp.', 'Spline interp.');

%% averaged result by method
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
plot(ms/n, sep_res);
legend('CCS', 'CS', 'Linear interp.', 'Spline interp.');
