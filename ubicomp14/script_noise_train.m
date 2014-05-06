clear
clc

s = load('d:\document\codes\data\NYC-noise-complaint\east_noise.txt');
dim = size(s);
idx = load('noise_top100.txt');

mode = 3; % 1 calc training set; 2 calc test set; 3 train

day1 = '12/31/2009';
day2 = '3/14/2014';
trainstop = '12/31/2012';
trainlength = daysact(day1, trainstop);
testlength = daysact(trainstop, day2);

intv = 365;
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
        days = daysact(trainstop,[int2str(s(i,4)) '/' int2str(s(i,5)) '/' int2str(s(i,3))]);
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
    end

elseif (mode == 3)
    b = zeros(100,100,15);
    data = load('noise_training_year.txt');
    for k = 1:15
        tb = orth(simple_ksvd(data, 100, k, 1));
        b(1:size(tb,1),1:size(tb,2),k) = tb;
        fprintf('%d is over, rank=%d\n', k, rank(tb));
    end
end