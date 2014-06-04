function y_gaussian = gaussian_disturb(y, theta);

n = length(y);
y_gaussian = zeros(size(y));
for i = 1:n
    y_gaussian(i) = normrnd(y(i), y(i) * theta);
end

    