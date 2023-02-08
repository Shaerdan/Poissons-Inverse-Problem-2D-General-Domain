clc; clear all; close all;
smoothsteps = 1;
rmax = 1;
rmin = 0.8;
n_mesh_bc = 400;
l_plot_smoothing = 1;
l_grad_test = 1;
curve_num = 1;
n_mesh_bulk = 400;
N_source = 4;
% N_ob = 20;
% N_skip = 5;
var_noise = 0;
% pattern = 3; % 1 for 'uniform partial', 2 for 'random', 3 for 'all'

% [H] = observationOperator(n_mesh_bc,N_ob,N_skip,pattern);

for itest = 1:curve_num
    [x_curve,y_curve] = randomsmoothcurveB2(smoothsteps,rmax,rmin,n_mesh_bc,l_plot_smoothing);
    hold on;
end
bc_curve = [x_curve;y_curve];
[x_rand,y_rand] = gen_mesh_points(x_curve,y_curve,n_mesh_bulk);
[x_source,y_source] = gen_mesh_points(x_curve,y_curve,N_source);

plot_num = 3;
figure(plot_num)
scatter(x_rand,y_rand,'.'); hold on; plot(x_curve,y_curve,'r*');
plot_num = plot_num +1;

[scaling_factor,w] = conformalscaling(x_curve,y_curve,n_mesh_bc);
figure(plot_num)
plot(scaling_factor,'DisplayName','Scaling factor')
xlabel('boundary mesh indx');
ylabel('w')
legend show
plot_num = plot_num + 1;

vertices = [x_curve, y_curve; x_rand, y_rand];
faces = convhull(x_curve, y_curve);
tri = delaunayTriangulation(vertices);
triplot(tri);hold on; scatter(x_source,y_source,'r*');
%K = stiffMatrix(tri);
% generate random intensity that sums to 0;
% S_true = randn(n_mesh_bulk,1); S_mean = mean(S_true);
% S_true = S_true - S_mean;
S_true = zeros(n_mesh_bulk,1);
S_true(1) = 1; S_true(floor(n_mesh_bulk/2)) = -1;
% sum(S_true)
Zi = [x_rand;y_rand];  z_ref = [x_curve;y_curve];
phi_bc_B2 =  soln_bc_B2(S_true,Zi,z_ref,n_mesh_bulk);
figure(plot_num)
plot(phi_bc_B2,'r-','Displayname','\phi_{B^{2}(0,1)}'); hold on;
phi_forward = scaling_factor.*phi_bc_B2;
plot(phi_forward,'k-','Displayname','\phi_{poly} = w.*\phi_{B^{2}(0,1)}');
xlabel('boundary mesh indx');
ylabel('\phi_{B^2(0,1)},\phi_{poly}');
legend show
legend show
plot_num = plot_num + 1;

figure(plot_num)
plot(phi_forward,'k-','Displayname','\phi_{poly}');
xlabel('boundary mesh indx');
ylabel('\phi_{poly}');
legend show
plot_num = plot_num + 1;

phi_data = phi_forward + var_noise*rand(n_mesh_bc,1);
figure(plot_num)
plot(phi_data,'k-','Displayname','\phi_{data}')
xlabel('boundary mesh indx');
ylabel('\phi_{poly}');
legend show
plot_num = plot_num + 1;




%gradient test:
if l_grad_test == 1
[f1,grad_f1] = costfun(S_true,Zi,z_ref,n_mesh_bulk,n_mesh_bc,scaling_factor,phi_data);
h = grad_f1/norm(grad_f1);
for itest = 1:14
    a(itest)= 10^(-itest);
    [f2,grad_f2] = costfun(S_true+a(itest)*h,Zi,z_ref,n_mesh_bulk,n_mesh_bc,scaling_factor,phi_data);
    PHI_a(itest) = abs(((f2-f1)/(a(itest)*h'*grad_f1))-1);
end
figure(plot_num)
loglog(a,PHI_a)
plot_num = plot_num + 1;
end

S0 = randn(n_mesh_bulk,1);
opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',5000,'MaxFunctionEvaluations',90000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-16,'CheckGradients',false,'SpecifyObjectiveGradient',true);
problem = createOptimProblem('fmincon','x0',S0, ...
    'objective',@(X) costfun(X,Zi,z_ref,n_mesh_bulk,n_mesh_bc,scaling_factor,phi_data),'options',opts);%...
    %'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,...
    %'options',opts);

S_soln = fmincon(problem);
plot_num = plot_num +1;
figure(plot_num)
bar(S_soln-S_true)
plot_num = plot_num +1;
figure(plot_num)
bar(S_soln); hold on; bar(S_true)
plot_num = plot_num +1;
n_source = 10;
X0 = randn(3*n_source,1);
opts = optimoptions(@fminunc,'PlotFcns',{@optimplotfval,@optimplotx},...
    'MaxIterations',5000,'MaxFunctionEvaluations',90000,...
    'OptimalityTolerance',1e-16,'StepTolerance',1e-16,'CheckGradients',true,'SpecifyObjectiveGradient',true);
problem = createOptimProblem('fmincon','x0',X0, ...
    'objective',@(X) costfun_nonlinear(X,z_ref,n_source,n_mesh_bc,scaling_factor,phi_data),'options',opts);%...
    %'Aeq',Aa,'beq',ba,'lb',lb,'ub',ub,...
    %'options',opts);
X_soln = fmincon(problem);
plot_num = plot_num +1;