mask = zeros(dim);
mask(samplex)=1;
HeatMap(mask);

% figure;
f = fopen('s_2d.txt', 'w');

figure;
hold on;

xidx = 0.05:0.05:0.8;

if (exist('res_ca', 'var'))
    [~, in] = sort(res_ca(:,1));
    res_ca = res_ca(in, :);
    s_ca = intervaled(res_ca, mfracn);
    plot(s_ca(:,6),s_ca(:,4), 'b-'); 
    [x y] = size(s_ca);
    fprintf(f, 's_ca\n');
    for i = 1 : x
        for j = 1 : y
            fprintf(f, '%f ', s_ca(i,j));
        end
        fprintf(f, '\n');
    end
end

if (exist('res_pair', 'var'))
    [~, in] = sort(res_pair(:,1));
    res_pair = res_pair(in, :);
    s_pair = intervaled(res_pair, mfracn);
    plot(s_pair(:,6),s_pair(:,4),'r--');
    [x y] = size(s_pair);
    fprintf(f, 's_pair\n');
    for i = 1 : x
        for j = 1 : y
            fprintf(f, '%f ', s_pair(i,j));
        end
        fprintf(f, '\n');
    end
end

if (exist('res_random', 'var')) 
    [~, in] = sort(res_random(:,1));
    res_random = res_random(in, :);
    s_random = intervaled(res_random, mfracn);
    plot(s_random(:,6),s_random(:,4),'m:'); 
    [x y] = size(s_random);
    fprintf(f, 's_random\n');
    for i = 1 : x
        for j = 1 : y
            fprintf(f, '%f ', s_random(i,j));
        end
        fprintf(f, '\n');
    end
end

if (exist('res_greedy', 'var')) 
    [~, in] = sort(res_greedy(:,1));
    res_greedy = res_greedy(in, :);
    s_greedy = intervaled(res_greedy, mfracn);
    plot(s_greedy(:,6),s_greedy(:,4),'c--.'); 
    [x y] = size(s_greedy);
    fprintf(f, 's_greedy\n');
    for i = 1 : x
        for j = 1 : y
            fprintf(f, '%f ', s_greedy(i,j));
        end
        fprintf(f, '\n');
    end
end

if (exist('res_local', 'var')) 
    [~, in] = sort(res_local(:,1));
    res_local = res_local(in, :);
    s_local = intervaled(res_local, mfracn);
    plot(s_local(:,6),s_local(:,4),'r--.'); 
    [x y] = size(s_local);
    fprintf(f, 's_greedy\n');
    for i = 1 : x
        for j = 1 : y
            fprintf(f, '%f ', s_local(i,j));
        end
        fprintf(f, '\n');
    end
end

fclose(f);