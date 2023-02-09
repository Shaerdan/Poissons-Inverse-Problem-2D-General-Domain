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
N_source = 4;

for itest = 1:curve_num
    [x_curve,y_curve] = randomsmoothcurveB2(smoothsteps,rmax,rmin,n_mesh_bc,l_plot_smoothing);
    hold on;
end
bc_curve = [x_curve;y_curve];
[x_source,y_source] = gen_mesh_points(x_curve,y_curve,N_source);

% generate random intensity that sums to 0;
S_true = randn(N_source,1); S_mean = mean(S_true);
S_true = S_true - S_mean;


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

tri = delaunayTriangulation([x_curve, y_curve]);

figure(901)
triplot(tri);hold on;
scatter(x_source,y_source,'r*');
var_noise = 0;
lambda = 0.1;
model = createpde();
stlwrite(tri,'tritext2D.stl','text')
geo = importGeometry(model,"tritext2D.stl");

g = @(region,state)gfun ;
mesh=generateMesh(model,"Hmax",.05,"GeometricOrder","linear");
figure(9090)
pdegplot(model.Geometry)
disp(strcat('Number of Nodes = ','',num2str(length(mesh.Nodes))));
applyBoundaryCondition(model,"neumann", ...
    "Face",1:model.Geometry.NumFaces, ...
    "g",g);
bndNodes0 = boundary(mesh.Nodes(1,:)',mesh.Nodes(2,:)');
bndNodes = unique(bndNodes0);


figure(1010)
scatter(mesh.Nodes(1,bndNodes)',mesh.Nodes(2,bndNodes)')

phi = @(xi, eta) [(1 - xi - eta), xi, eta];
conv_kernel = [];
tol = 1e-5;
i_sigma = 1;
eps = 1*1e-1;
f = @(location,state)fsource(location,state,S_true,x_source,y_source,N_source,eps,conv_kernel);
specifyCoefficients(model,"m",0,...
    "d",0,...
    "c",1,...
    "a",0,...
    "f",f);
FEMn = assembleFEMatrices(model,"nullspace");
noise = sqrt(var_noise)*randn(size(bndNodes));

cond(FEMn.Kc)

tikh = FEMn.Kc'*FEMn.Kc + lambda*speye(size(FEMn.Kc));
FcKinv = tikh\FEMn.Kc'*FEMn.Fc;
soln_FE_n_tikh(:,1) = FEMn.B*(FcKinv) + FEMn.ud;
Data_Bc = soln_FE_n_tikh(bndNodes,1) + noise;
figure(10)
plot(Data_Bc)

if l_plot_loop == 1
    figure(plot_num)
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
X0 = rand(size(soln_FE_n_tikh));

options0 = optimoptions('fmincon','CheckGradients',false,'SpecifyObjectiveGradient',false,...
    'PlotFcn','optimplotfval','MaxIterations',20000,'MaxFunctionEvaluations',1e5);
soln_reconstructed = fmincon(@(X)costfun_FEM(X,A,Data_Bc,bndNodes,lambda),X0,[],[],[],[],[],[],[],options0);
figure(651)
Data_Estimated = A*soln_reconstructed;
plot(Data_Estimated(bndNodes))

figure(602)
scatter(model.Geometry.Vertices(:,1),model.Geometry.Vertices(:,2))

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
