function [phi_bc] = soln_bc_B2(S,Zi,z_ref,N_source)
% function to compute the forward solution in cartesian coordinates
% Zi is the source location, z_ref is the boundary reference point
% locations.

n_ref = length(z_ref)/2;

x_ref = z_ref(1:n_ref)'; y_ref = z_ref(n_ref+1:end)';
Xi = Zi(1:N_source); Yi = Zi(N_source + 1:end);
phi = zeros(n_ref,N_source);
for i = 1:N_source
    for j=1:n_ref
        dx_reftoSource = Xi(i) - x_ref(j);
        dy_reftoSource = Yi(i) - y_ref(j);
        phi(j,i) = (S(i)/(2*pi))*log(norm(dx_reftoSource)+ norm(dy_reftoSource));
    end
end
phi_bc = 2*sum(phi,2);
end

