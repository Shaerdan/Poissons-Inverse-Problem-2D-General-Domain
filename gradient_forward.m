function [g] = gradient_forward(scaling_factor,phi,N_source,n_ref)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
g = zeros(n_ref,N_source);
for i = 1:N_source
g(:,i) = scaling_factor.*phi(:,i);
end
g = 2*g;
end

