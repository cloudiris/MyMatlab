clear
clc

% s = load('d:\document\codes\data\NYC-graffiti\east_graffiti.txt');
s = load('d:\document\codes\data\NYC-rat\east_rat.txt');
dim = size(s);
idx = load('top100.txt');

mode = 1; % 1 calc training set; 2 calc test set; 3 train

intv = 365;
day1 = '12/31/2009';
day2 = '3/12/2014';
stoptrain = '12/31/2012';
trainlength = daysact(day1, stoptrain);
testlength = daysact(stoptrain, day2);

training = zeros(100, trainlength-intv+1);
test = zeros(100, testlength-intv+1);
t1 = 0; t2 = 0; t3 = 0;

if (mode == 1)
    for i = 1 : length(s)
        if (s(i,3) >= 2013) continue; end
        tic;
        days = daysact(day1,[int2str(s(i,4)) '/' int2str(s(i,5)) '/' int2str(s(i,3))]);
        t1 = t1 + toc;
        tic;
        intv_start = max(1, days-intv+1);
        intv_end = min(trainlength-intv+1, days);
        intvs = intv_start : intv_end;
        x = ceil((s(i,1) - 970000)/5000);
        y = ceil((s(i,2) - 140000)/7000);
        coor = (y - 1)*20 + x;
        offset = find(idx==coor);
        t2 = t2 + toc;
        tic;
        if (~isempty(offset))
            intvn = length(intvs);
            for intvi = 1:intvn
                in = intvs(intvi);
                training(offset, in) = training(offset, in) + 1;
            end
        end
        t3 = t3 + toc;
    end
    
elseif (mode == 2)
    for i = 1 : length(s)
        if (s(i,3) < 2013) continue; end
        tic;
        days = daysact(stoptrain,[int2str(s(i,4)) '/' int2str(s(i,5)) '/' int2str(s(i,3))]);
        t1 = t1 + toc;
        tic;
        intv_start = max(1, days-intv+1);
        intv_end = min(testlength-intv+1, days);
        intvs = intv_start : intv_end;
        x = ceil((s(i,1) - 970000)/5000);
        y = ceil((s(i,2) - 140000)/7000);
        coor = (y - 1)*20 + x;
        offset = find(idx==coor);
        t2 = t2 + toc;
        tic;
        if (~isempty(offset))
            intvn = length(intvs);
            for intvi = 1:intvn
                in = intvs(intvi);
                test(offset, in) = test(offset, in) + 1;
            end
        end
        t3 = t3 + toc;
        fprintf('%d\n',i);
    end

elseif (mode == 3)
    kk = 15;
    b = zeros(100,100,kk);
%     data = load('rat_training_100_209.txt');
    data = load('rat_training_100_year.txt');
    err = zeros(kk, 1);
    for k = 1:kk
        [tb, o] = simple_ksvd(data, 100, k, 1);
        tb = orth(tb);
        err(k) = o.besterr;
        b(1:size(tb,1),1:size(tb,2),k) = tb;
        fprintf('%d is over\n', k);
    end
end