function [x_curve,y_curve] = randomsmoothcurveB2(smoothsteps,rmax,rmin,n_mesh,l_plot_smoothing)
%   This function generate a random 2D curve and smooth it
%   x_curve: x coordinates of the curve
%   y_curve: y coordinates of the curve
%   radius is generated from uniform random distribution of the interval
%   [rmin rmax]

radius_thetai = (rmax-rmin)*rand(n_mesh,1) + rmin;
radius_thetai(end) =  radius_thetai(1); 
theta_mesh = linspace(0,2*pi,n_mesh)';
y_curve = radius_thetai.*sin(theta_mesh);
x_curve = radius_thetai.*cos(theta_mesh);
if l_plot_smoothing == 1
figure(1)
plot(x_curve,y_curve)
end
for i = 1:smoothsteps
x_curve = smooth(x_curve,'sgolay');
y_curve = smooth(y_curve,'sgolay');
x_curve(end) = x_curve(1); y_curve(end) = y_curve(1);
end
if l_plot_smoothing == 1
figure(2)
plot(x_curve,y_curve);
end
end

