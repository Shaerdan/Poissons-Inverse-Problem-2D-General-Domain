clc; clear all; close all;
smoothsteps = 1;
rmax = 1;
rmin = 0.8;
n_mesh_bc = 180;
l_plot_smoothing = 1;
l_plot_conformal = 0;
l_plot_loop = 0;
l_grad_test = 1;
curve_num = 1;
n_mesh_bulk = 400;
N_source = 10;
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
% [x_rand,y_rand] = gen_mesh_points(x_curve,y_curve,n_mesh_bulk);
[x_source,y_source] = gen_mesh_points(x_curve,y_curve,N_source);

% plot_num = 3;
% figure(plot_num)
% scatter(x_rand,y_rand,'.'); hold on; plot(x_curve,y_curve,'r*');
% plot_num = plot_num +1;

% generate random intensity that sums to 0;
S_true = randn(N_source,1); S_mean = mean(S_true);
S_true = S_true - S_mean;
% S_true = zeros(n_mesh_bulk,1);
% S_true(1) = 1; S_true(floor(n_mesh_bulk/2)) = -1;
% sum(S_true)
% Zi = [x_rand;y_rand];  z_ref = [x_curve;y_curve];


if l_plot_conformal == 1
    [scaling_factor,w] = conformalscaling(x_curve,y_curve,n_mesh_bc);
    
    phi_bc_B2 =  soln_bc_B2(S_true,Zi,z_ref,N_source);
    phi_forward = scaling_factor.*phi_bc_B2;
    figure(plot_num)
    plot(phi_forward,'k-','Displayname','\phi_{poly}');
    xlabel('boundary mesh indx');
    ylabel('\phi_{poly}');
    legend show
    plot_num = plot_num + 1;
    
end
% vertices = [x_curve, y_curve; x_rand, y_rand];
% faces = convhull(x_curve, y_curve);
% tri = delaunayTriangulation(vertices);
tri = delaunayTriangulation([x_curve, y_curve]);
% filename1 = 'bc_geometry';
% stlwrite(filename1, tri)% tri = [x_curve, y_curve];
figure(901)
triplot(tri);hold on;
scatter(x_source,y_source,'r*');


model = createpde();
stlwrite(tri,'tritext2D.stl','text')
geo = importGeometry(model,"tritext2D.stl");

% g = geometryFromEdges(model,@circleg);
g = @(region,state)gfun ;
mesh=generateMesh(model,"Hmax",.05,"GeometricOrder","linear");

disp(strcat('Number of Nodes = ','',num2str(length(mesh.Nodes))));
applyBoundaryCondition(model,"neumann", ...
    "Face",1:model.Geometry.NumFaces, ...
    "g",g);
Nf2 = 1:64;
elemIDs = 1:60;
figure(1010)
scatter(mesh.Nodes(1,elemIDs)',mesh.Nodes(2,elemIDs)')
% phi = @(xi, eta) [1-xi, 1+xi, 1-eta, 1+eta];
% phi = @(xi, eta) [(1 - xi) * (1 - eta) / 4, (1 + xi) * (1 - eta) / 4,...
%     (1 + xi) * (1 + eta) / 4, (1 - xi) * (1 + eta) / 4];
phi = @(xi, eta) [(1 - xi - eta), xi, eta];
% gaussKernel = fspecial('gaussian', [5 5], 1);
conv_kernel = [];
lambda = 0.1;
tol = 1e-5;
i_sigma = 1;
eps = 5*1e-2;
f = @(location,state)fsource(location,state,S_true,x_source,y_source,N_source,eps,conv_kernel);
specifyCoefficients(model,"m",0,...
    "d",0,...
    "c",1,...
    "a",0,...
    "f",f);
FEMn = assembleFEMatrices(model,"nullspace");
cond(FEMn.Kc)
%     FEMs = assembleFEMatrices(model,"stiff-spring");
%     soln_FE_n = FEMn.B*(FEMn.Kc\FEMn.Fc) + FEMn.ud;
%     Kinv = ridge(FEMn.Kc, lambda);
tikh = FEMn.Kc'*FEMn.Kc + lambda*speye(size(FEMn.Kc));
FcKinv = tikh\FEMn.Kc'*FEMn.Fc;
soln_FE_n_tikh(:,1) = FEMn.B*(FcKinv) + FEMn.ud;
%     noise = randn(size(FEMn.Fc));
noise = 0;
Data_Fc = FEMn.Fc + noise;
%     soln_FE_s = FEMs.Ks\FEMs.Fs;
results = solvepde(model,phi);

% extract boundary solution:
phi_forward_FE_raw(:,1) = results.NodalSolution(Nf2);
phi_forward_FE_tikh(:,1) = soln_FE_n_tikh(Nf2,1);

soln_FE_raw(:,1) = results.NodalSolution;
if l_plot_loop == 1
    figure(plot_num)
    pdeplot(model,"XYData",results.NodalSolution)
    legend('unregularised')
    legend show
    plot_num = plot_num + 1;
    
    figure(plot_num)
    pdeplot(model,"XYData",soln_FE_n_tikh)
    legend('Tikhonov on stiffness matrix')
    legend show
    
    plot_num = plot_num + 1;
    
    figure(600)
    plot(phi_forward_FE_raw); hold on;
end

A = FEMn.Kc;
% A_Tikh = FEMn.Kc'*FEMn.Kc + lambda*speye(size(FEMn.Kc));
% Data_Fc_Tikh = FEMn.Kc'*Data_Fc;
X0 = rand(size(Data_Fc));
figure(650)
lambda = 0.1;
plot(Data_Fc)
options0 = optimoptions('fmincon','CheckGradients',false,'SpecifyObjectiveGradient',false,...
    'PlotFcn','optimplotfval','MaxIterations',20000,'MaxFunctionEvaluations',1e5);
soln_reconstructed = fmincon(@(X)costfun_FEM(X,A,Data_Fc,lambda),X0,[],[],[],[],[],[],[],options0);
figure(651)
plot(A*soln_reconstructed)
% soln_reconstructed_Tikh = fmincon(@(X)costfun_FEM(X,A_Tikh,Data_Fc_Tikh),X0,[],[],[],[],[],[],[],options0);
% figure(1021)
% pdeplot(model,"XYData",soln_reconstructed_Tikh(:,1),"FaceAlpha",0.9); hold on;

figure(600)
plot(phi_forward_FE_raw(:,1),'DisplayName',strcat('Unregularised, eps =',eps)); hold on;

xlabel('boundary nodes')
ylabel('\phi(\partial \Omega)')
legend show


figure(601)
plot(phi_forward_FE_tikh(:,1),'DisplayName',strcat('Tikhonov regularised, eps =',eps)); hold on;
xlabel('boundary nodes')
ylabel('\phi(\partial \Omega)')
legend show

figure(602)
scatter(model.Geometry.Vertices(:,1),model.Geometry.Vertices(:,2))

figure(603)
pdeplot(model,"XYData",soln_FE_raw(:,1))
xlabel('x')
ylabel('y')
legend('unregularised')
legend show

figure(604)
markersize_scaleup = 100;
pdeplot(model,"XYData",soln_FE_n_tikh(:,1),"FaceAlpha",0.9); hold on;
S1 = scatter(x_source(S_true<=0),y_source(S_true<=0),markersize_scaleup*abs(S_true(S_true<=0)),'ko'); hold on;
S2 = scatter(x_source(S_true>0),y_source(S_true>0),markersize_scaleup*abs(S_true(S_true>0)),'k+');
xlabel('x')
ylabel('y')
legend('Tikhonov regularised')
legend show

figure(1020)
pdeplot(model,"XYData",soln_reconstructed(:,1),"FaceAlpha",0.9); hold on;
xlabel('x')
ylabel('y')
legend('Inverse Solution')
legend show
