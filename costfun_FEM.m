function [f,g] = costfun_FEM(X,A,Data,bndNodes,lambda)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
u = A*X;
u_traced = u(bndNodes);
f = 0.5d0*norm(u_traced-Data)^2 + 0.5d0*lambda*norm(X,2)^2; %+ (lambda/100)*norm(X,1);
% g_tikh = (A')*(A*X-Data) + lambda*X;

g= [];
end

