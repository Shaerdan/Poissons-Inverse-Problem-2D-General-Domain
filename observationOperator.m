function [H] = observationOperator(N_bc,N_ob,N_skip,pattern)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
H = zeros(N_bc,1);
if pattern == 1
    rand_int = randperm(N_bc);
    index = rand_int(1:N_ob);
elseif pattern == 2
    index = 1:N_skip:N_bc;
elseif pattern == 3
    index = 1:N_bc;
end
for i = index
    H(i) = 1;
end

end

