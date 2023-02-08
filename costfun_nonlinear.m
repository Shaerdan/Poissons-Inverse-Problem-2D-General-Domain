function [f,grad_f] = costfun_nonlinear(X,z_ref,N_source,n_mesh_bc,scaling_factor,phi_data)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
S = X(1:N_source);
Zi = X(N_source+1:end);
[phi_bc_B2,phi_g,phi_g_loc_x,phi_g_loc_y] =  soln_bc_B2(S,Zi,z_ref,N_source);
phi_forward = scaling_factor.*phi_bc_B2;
descrep = phi_forward-phi_data;
f = 0.5d0*norm(descrep)^2;
g_S = gradient_forward(scaling_factor,phi_g,N_source,n_mesh_bc);
g_loc = gradient_forward_loc(scaling_factor,phi_g_loc_x,phi_g_loc_y,N_source,n_mesh_bc);
g = [g_S g_loc];
grad_f = (g'*descrep);
% grad_f = [];
end

