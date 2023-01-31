clc; clear all; close all;
smoothsteps = 20;
rmax = 1;
rmin = 0.7;
n_mesh_bc = 400;
l_plot_smoothing = 1;
curve_num = 1;
n_mesh_bulk = 40;
for itest = 1:curve_num
[x_curve,y_curve] = randomsmoothcurveB2(smoothsteps,rmax,rmin,n_mesh_bc,l_plot_smoothing);
hold on;
end
bc_curve = [x_curve;y_curve];
[x_rand,y_rand] = gen_mesh_points(x_curve,y_curve,n_mesh_bulk);
figure(3)
scatter(x_rand,y_rand,'.'); hold on; plot(x_curve,y_curve,'r*');

[scaling_factor] = conformalscaling(x_curve,y_curve,n_mesh_bc);
