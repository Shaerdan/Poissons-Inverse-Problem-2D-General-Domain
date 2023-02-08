function [g] = gradient_forward_loc(scaling_factor,phi_x,phi_y,N_source,n_ref)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
g = zeros(n_ref,2*N_source);
for i = 1:N_source
g(:,i) = scaling_factor.*phi_x(:,i);
end
for i = N_source+1:2*N_source
g(:,i) = scaling_factor.*phi_y(:,i-N_source);
end
g = 2*g;
end

