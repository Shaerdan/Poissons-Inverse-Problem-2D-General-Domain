function [f,g] = costfun_FEM(X,A,Data,bndNodes,lambda,reg_method)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
u_bc = A*X;

% if reg_method == 1
%     f = 0.5d0*norm(u_traced-Data)^2 + 0.5d0*lambda*norm(X,2)^2; %+ (lambda/100)*norm(X,1);
% else
    f = 0.5d0*norm(u_bc-Data)^2 + 0.5d0*lambda*norm(X,1);
% end
% g_tikh = (A')*(A*X-Data) + lambda*X;

g= [];
end

