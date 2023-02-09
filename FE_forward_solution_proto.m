clc; clear all; close all;
smoothsteps = 1;
rmax = 1;
rmin = 0.7;
n_vertices_bc_geometry = 380;
H_max_FEM = 0.5;
eps = 5*1e-2; % margin for the Kronecker delta function
l_plot_smoothing = 1;
l_plot_conformal = 0;
l_plot_loop = 0;
l_grad_test = 1;
curve_num = 1;
n_mesh_bulk = 400;
N_source = 2;
var_noise = 0;
lambda_FEM = 0.1;
lambda_inv_tikh = 5*1e-1;
lambda_inv_lasso = 1e-5;
reg_method = 2; % 1 for Tikhonov, 2 for lasso
if reg_method == 1
    lambda_inv = lambda_inv_tikh;
elseif reg_method == 2
    lambda_inv = lambda_inv_lasso;
end

for itest = 1:curve_num
    [x_curve,y_curve] = randomsmoothcurveB2(smoothsteps,rmax,rmin,n_vertices_bc_geometry,l_plot_smoothing);
    hold on;
end
bc_curve = [x_curve;y_curve];
[x_source,y_source] = gen_mesh_points(x_curve,y_curve,N_source);

% generate random intensity that sums to 0;
S_true = randn(N_source,1); S_mean = mean(S_true);
S_true = S_true - S_mean;


if l_plot_conformal == 1
    [scaling_factor,w] = conformalscaling(x_curve,y_curve,n_vertices_bc_geometry);
    
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

model = createpde();
stlwrite(tri,'tritext2D.stl','text')
geo = importGeometry(model,"tritext2D.stl");

g = @(region,state)gfun ;
% mesh=generateMesh(model,"Hmax",H_max_FEM,"GeometricOrder","linear");
mesh=generateMesh(model,"Hgrad",1.1,"GeometricOrder","linear");

figure(9090)
[p,e,t] = meshToPet(mesh);
pdeplot(p,e,t)
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
f = @(location,state)fsource(location,state,S_true,x_source,y_source,N_source,eps,conv_kernel);
specifyCoefficients(model,"m",0,...
    "d",0,...
    "c",1,...
    "a",0,...
    "f",f);
FEMn = assembleFEMatrices(model,"nullspace");
noise = sqrt(var_noise)*randn(size(bndNodes));

cond(FEMn.Kc)

tikh = FEMn.Kc'*FEMn.Kc + lambda_FEM*speye(size(FEMn.Kc));
FcKinv = tikh\FEMn.Kc'*FEMn.Fc;
soln_FE_n_tikh(:,1) = FEMn.B*(FcKinv) + FEMn.ud;

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
A_tikh_inv = tikh\FEMn.Kc';
nndoes = length(A_tikh_inv); nbcnodes = length(bndNodes);
H_obs = zeros(nbcnodes,nndoes);
H_obs_diag = ones(nbcnodes,1); % observation pattern on the boundary
H_obs(1:nbcnodes,1:nbcnodes) = diag(H_obs_diag);
H_obs = sparse(H_obs);
A_forward = H_obs*A_tikh_inv;
X0 = 1e-3*ones(size(soln_FE_n_tikh));
Data_Bc = H_obs*soln_FE_n_tikh + noise;
figure(10)
plot(Data_Bc)
options0 = optimoptions('fmincon','CheckGradients',false,'SpecifyObjectiveGradient',false,...
    'PlotFcn','optimplotfval','MaxIterations',1e7,'MaxFunctionEvaluations',1e7);
sourceterm_reconstructed = fmincon(@(X)costfun_FEM(X,A_forward,Data_Bc,bndNodes,lambda_inv,reg_method),X0,[],[],[],[],[],[],[],options0);
figure(651)
Soln_Inv = A_tikh_inv*sourceterm_reconstructed;
plot(Soln_Inv(bndNodes))

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
pdeplot(model,"XYData",Soln_Inv(:,1),"FaceAlpha",0.9); hold on;
xlabel('x')
ylabel('y')
legend('Inverse Solution')
legend show

figure(1021)
pdeplot(model,"XYData",sourceterm_reconstructed(:,1),"FaceAlpha",0.9); hold on;
xlabel('x')
ylabel('y')
legend('Inverse Solution')
legend show
