function [f,grad_f] = costfun(S,Zi,z_ref,n_mesh_bulk,n_mesh_bc,scaling_factor,phi_data)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

[phi_bc_B2,phi_g] =  soln_bc_B2(S,Zi,z_ref,n_mesh_bulk);
phi_forward = scaling_factor.*phi_bc_B2;
descrep = phi_forward-phi_data;
f = 0.5d0*norm(descrep)^2;
[g] = gradient_forward(scaling_factor,phi_g,n_mesh_bulk,n_mesh_bc);
grad_f = (g'*descrep);
% grad_f = [];
end

