b = load('battery_level_w2.txt');
ib = 1 - b;

h = 10000;
l = 1;

li = (h-l) * ib + l;
ex = l*(h/l).^ib;

di = 50;
lo = (h-l)*log(ib*(di - 1) + 1)/log(di) + l;

save('battercost_li_w2(big).txt','li','-ascii');
save('battercost_ex_w2(big).txt','ex','-ascii');
save('battercost_lo_w2(big).txt','lo','-ascii');

figure;

[cc,b] = hist(li(:),100);
cc = cc / numel(li);
subplot(3,1,1);
bar(b,cc);
set(gca,'linewidth',2);

[cc,b] = hist(ex(:),100);
cc = cc / numel(li);
subplot(3,1,2);
bar(b,cc);
set(gca,'linewidth',2);

[cc,b] = hist(lo(:),100);
cc = cc / numel(li);
subplot(3,1,3);
bar(b,cc);
set(gca,'linewidth',2);