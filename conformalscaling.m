function [scaling_factor,w] = conformalscaling(x,y,N)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% Define the vertices of the polygon


% Create a set of points in the unit disk
t = linspace(0, 2*pi, N);
z = exp(1i * t);
% [X, Y] = meshgrid(linspace(-1, 1, N), linspace(-1, 1, N));
% z = X + 1i * Y;
% idx = abs(z) <= 1;
% z = z(idx);
% Compute the conformal mapping using the Schwarz-Christoffel formula
w = boundary_integral_equation(x, y, z);
% w = lft(x, y, z);
% Calculate the scaling factor at each point in the disk
scaling_factor = zeros(N, 1);
for i = 1:N
    dw_dz = (w(mod(i, N) + 1) - w(i)) / (z(mod(i, N) + 1) - z(i));
    scaling_factor(i) = abs(dw_dz / (z(i)));
end
end

