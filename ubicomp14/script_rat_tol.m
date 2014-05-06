%% tol ~ acc

tols = 0.1:0.05:0.5;
toln = length(tols);
acc = zeros(mn, toln, methodn);
    
sep_res = zeros(mn, methodn);
for mi = 1 : mn
    for methodi = 1 : methodn
        for toli = 1 : toln
            tol = tols(toli);
        
            sep_res(mi, methodi) = 0;
            for si = 1 : size(test,2)
                s = test(:,si);
                tolv = s.*tol;
    %             sep_res(mi, methodi) = sep_res(mi, methodi) + snr(s, re_by_method(:,si,mi,methodi));
                sep_res(mi, methodi) = sep_res(mi, methodi) + sum(abs(s-re_by_method(:,si,mi,methodi))<tolv)/n;
            end
            sep_res(mi, methodi) = sep_res(mi, methodi) / size(test,2);
            acc(mi, toli, methodi) = sep_res(mi, methodi);
        end
    end
end