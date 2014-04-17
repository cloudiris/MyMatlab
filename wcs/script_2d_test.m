% c = exprnd(11.9, 64,64) + 1;
% c(find(c>100)) = 100;
c = positioncost;
cv = hilbertcurve(c);
sc = sort(cv);
mask = ones(1, 4096);
c2v = zeros(1, 4096);
c2 = zeros(64);

d = 110;
target = 1190;

for i = 1 : 4096
    [~, idx] = min(abs(sc-target));
    c2v(i) = sc(idx);
    target = abs(normrnd(sc(idx), d));
    sc(idx) = inf;
end

% c2 = izigzag(c2v, 64, 64);
c2 = antihc(c2v);
surf(c2);
% plot(c2v);