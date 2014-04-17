file = ['cc0518.txt';'cc0519.txt';'cc0520.txt';'cc0521.txt';'cc0522.txt'];
ccount = [];
for i = 1 : 5
    filename = ['taxidata\' file(i,:)];
    cc = load(filename);
    for j = 1 : length(cc);
        if (length(ccount) < cc(j,1))
            ccount(cc(j,1)) = cc(j,2);
        else
            ccount(cc(j,1)) = ccount(cc(j,1)) + cc(j,2);
        end
    end
end

f = fopen('cabcount_w1.txt', 'w');
for i = 1 : numel(ccount)
    fprintf(f, '%d\n', ccount(i));
end
fclose(f);