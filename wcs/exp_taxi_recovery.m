clear
clc

workpath = 'd:\document\codes\matlab\wcs\';
seg = load('taxidata\segments.txt');
gt_all = load([workpath 'week_gt.txt']);
m_all = load([workpath 'week_m_ca_0.2.txt']);
mask_all = load([workpath 'week_mask_ca_0.2.txt']);

si1 = 401;
si2 = 428;

n = length(seg);
intvn = 26 * 7;

rank = 70;
t = 500;
lamda = 100;

disp('begin recovery\n');

tic;
    %% matrix completion
    r = matrixcompletion(m_all, mask_all, rank, lamda, t);

    %% inexact ALM MC
    
    %% time series
    len = intvn * 20;
    base = base_dct4(len);
    samplex = [];
    y = [];
    k = floor(0.01 * len);
    count = 0;
    for i = 1 : 20
        for j = 1 : intvn
            if (mask_all(j, i) == 1)
              count = count + 1;
              samplex(count) = (i - 1) * intvn + j;
              y(count) = m_all(j, i);
            end
        end
    end
    r = my_csf(y', samplex, base, k, 'OMP', 1);
toc;