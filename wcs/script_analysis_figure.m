figure;
plot(s_r(:,1),s_r(:,3), 'm:'); l1='Random';
hold on
plot(s_ca1(:,1),s_ca1(:,3), 'b^-'); l2 = 'CA(1)';
plot(s_cadot1(:,1),s_cadot1(:,3), 'b+-'); l3 = 'CA(0.1)';
plot(s_ca5(:,1),s_ca5(:,3), 'bd-'); l4 = 'CA(5)';

plot(s_2group(:,1),s_2group(:,3), 'r^--'); l5 = '2-Group';
plot(s_4group(:,1),s_4group(:,3), 'rd--'); l6 = '4-Group';

xlim([0 1]);
ylim([0 1]);


legend(l1, l2, l3, l4, l5, l6);