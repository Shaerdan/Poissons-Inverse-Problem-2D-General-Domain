function [f,g] = costfun_FEM(X,A,Data_Fc,lambda)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
f = 0.5d0*norm(A*X-Data_Fc)^2 + 0.5d0*lambda*norm(X,2)^2 %+ (lambda/100)*norm(X,1);
g = [];
end

