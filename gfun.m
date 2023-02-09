function [g] = gfun(region,state)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
g = region.^2 + state;
disp('region',region);
end

