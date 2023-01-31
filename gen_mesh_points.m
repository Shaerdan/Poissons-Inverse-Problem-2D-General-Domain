function [x_rand,y_rand] = gen_mesh_points(x_curve,y_curve,n_points)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% Define the number of random points to generate


% Define the bounding box that encloses the curve
xmin = min(x_curve);
xmax = max(x_curve);
ymin = min(y_curve);
ymax = max(y_curve);

% Initialize arrays to store the random points
x_rand = zeros(n_points, 1);
y_rand = zeros(n_points, 1);

% Set a counter to keep track of the number of points generated inside the curve
count = 0;

% Use a while loop to generate random points until the desired number of points inside the curve is reached
while count < n_points
    % Generate a random point within the bounding box
    x_rand_temp = xmin + (xmax-xmin)*rand();
    y_rand_temp = ymin + (ymax-ymin)*rand();
    
    % Check if the point is inside the curve using the ray casting method
    line = [x_rand_temp, y_rand_temp; x_rand_temp + 1e6, y_rand_temp];
    [xi, yi] = polyxpoly(x_curve, y_curve, line(:,1), line(:,2));
    if mod(length(xi), 2) == 1
        % If the point is inside the curve, increment the counter and store the point
        count = count + 1;
        x_rand(count) = x_rand_temp;
        y_rand(count) = y_rand_temp;
    end
end
end

