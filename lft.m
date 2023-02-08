function w = lft(x, y, z)
% LFT Computes the conformal mapping from the unit disk to a polygon using linear fractional transformation
% x, y: vertices of the polygon
% z: points in the unit disk
% w: image of the points in the disk in the complex plane

% Number of vertices of the polygon
n = length(x);

% Number of points in the disk
m = length(z);

% Allocate memory for the image of the points in the disk
w = zeros(m, 1);

% Compute the coefficients of the LFT
p = poly(z);
a = polyval(p, x + 1i * y);
b = polyval(p, 1 ./ conj(x + 1i * y));
fprintf('Number of elements in a: %d\n', length(a));
fprintf('Number of elements in b: %d\n', length(b));
fprintf('Number of elements in z: %d\n', length(z));

% Compute the LFT for each point in the disk
for i = 1:m
w(i) = (sum(a.*z(i)) + sum(b)) / (sum(conj(a).*z(i)) + sum(conj(b)));
end
end