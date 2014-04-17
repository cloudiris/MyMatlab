rn = zeros(130, 1);
for i = 1 : 130
    file = sprintf('taxidata\\w2_intv%d.txt',i);
    nrec = load(file);
    rn(i) = length(nrec);
    fprintf('intv%d: %d recs\n', i, rn(i));
end