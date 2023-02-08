function [phi_bc,phi_g,phi_g_loc_x,phi_g_loc_y] = soln_bc_B2(S,Zi,z_ref,N_source)
% function to compute the forward solution in cartesian coordinates
% Zi is the source location, z_ref is the boundary reference point
% locations. Also output phi_g for gradient computation.

n_ref = length(z_ref)/2;

x_ref = z_ref(1:n_ref)'; y_ref = z_ref(n_ref+1:end)';
Xi = Zi(1:N_source); Yi = Zi(N_source + 1:end);
phi = zeros(n_ref,N_source);
phi_g = phi;
phi_g_loc_x= phi;
phi_g_loc_y = phi;
for i = 1:N_source
    for j=1:n_ref
        dx_reftoSource = Xi(i) - x_ref(j);
        dy_reftoSource = Yi(i) - y_ref(j);
        phi_g(j,i) = (1/(2*pi))*log(dx_reftoSource^2+ dy_reftoSource^2);
        phi_g_loc_x(j,i) = (S(i)/(2*pi))*(1/(dx_reftoSource^2+ dy_reftoSource^2))*...
                            (2*dx_reftoSource);
        phi_g_loc_y(j,i) = (S(i)/(2*pi))*(1/(dx_reftoSource^2+ dy_reftoSource^2))*...
                            (2*dy_reftoSource);                        
        phi(j,i) = (S(i)*phi_g(j,i));
    end
end
phi_bc = 2*sum(phi,2);
end

