function [f] = fsource(location,state,S,x_source,y_source,nsource,eps,kernel)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% distance = zeros(nsource,1);
for i= 1:nsource
    distance = sqrt((location.x-x_source(i)).^2 + (location.y-y_source(i)).^2);
    %distance(distance<=0.01) = 0;
    f_i(:,i) = S(i)*exp(-(distance.^2)./(2*eps^2));
end
f = sum(f_i,2);
f = f';

if kernel ~= []
    % Reshape the source vector into a 2D matrix
    f_matrix = reshape(f, [length(f), 1]);
    
    % Perform the 2D convolution on the reshaped matrix
    f_conv = conv2(f_matrix, kernel, 'same');
    f = f_conv';
    % f_conv = reshape(f_matrix_conv, [length(f), 1]);
end
% f = location.x + location.y;
end

