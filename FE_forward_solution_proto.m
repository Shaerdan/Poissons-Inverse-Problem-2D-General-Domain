clc; clear all; close all;
smoothsteps = 1;
rmax = 1;
rmin = 0.8;
n_mesh_bc = 40;
l_plot_smoothing = 1;
l_plot_conformal = 1;
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
[x_rand,y_rand] = gen_mesh_points(x_curve,y_curve,n_mesh_bulk);
[x_source,y_source] = gen_mesh_points(x_curve,y_curve,N_source);

plot_num = 3;
figure(plot_num)
scatter(x_rand,y_rand,'.'); hold on; plot(x_curve,y_curve,'r*');
plot_num = plot_num +1;

% generate random intensity that sums to 0;
S_true = randn(N_source,1); S_mean = mean(S_true);
S_true = S_true - S_mean;
% S_true = zeros(n_mesh_bulk,1);
% S_true(1) = 1; S_true(floor(n_mesh_bulk/2)) = -1;
% sum(S_true)
Zi = [x_rand;y_rand];  z_ref = [x_curve;y_curve];


if l_plot_conformal == 1
    [scaling_factor,w] = conformalscaling(x_curve,y_curve,n_mesh_bc);
%     figure(plot_num)
%     plot(scaling_factor,'DisplayName','Scaling factor')
%     xlabel('boundary mesh indx');
%     ylabel('w')
%     legend show
%     plot_num = plot_num + 1;
    
    phi_bc_B2 =  soln_bc_B2(S_true,Zi,z_ref,N_source);
    phi_forward = scaling_factor.*phi_bc_B2;
    figure(plot_num)
    plot(phi_forward,'k-','Displayname','\phi_{poly}');
    xlabel('boundary mesh indx');
    ylabel('\phi_{poly}');
    legend show
    plot_num = plot_num + 1;    
%     figure(plot_num)
%     plot(phi_bc_B2,'r-','Displayname','\phi_{B^{2}(0,1)}'); hold on;
%     plot(phi_forward,'k-','Displayname','\phi_{poly} = w.*\phi_{B^{2}(0,1)}');
%     xlabel('boundary mesh indx');
%     ylabel('\phi_{B^2(0,1)},\phi_{poly}');
%     legend show
%     legend show
%     plot_num = plot_num + 1;
    
%     phi_data = phi_forward + var_noise*rand(n_mesh_bc,1);
%     figure(plot_num)
%     plot(phi_data,'k-','Displayname','\phi_{data}')
%     xlabel('boundary mesh indx');
%     ylabel('\phi_{poly}');
%     legend show
%     plot_num = plot_num + 1;
    
end
vertices = [x_curve, y_curve; x_rand, y_rand];
faces = convhull(x_curve, y_curve);
% tri = delaunayTriangulation(vertices);
tri = delaunayTriangulation([x_curve, y_curve]);
% tri = [x_curve, y_curve];
figure(901)
triplot(tri);hold on;
scatter(x_source,y_source,'r*');
%K = stiffMatrix(tri);



model = createpde();
stlwrite(tri,'tritext2D.stl','text')
importGeometry(model,"tritext2D.stl");
gfun = @(region,state)0;
mesh=generateMesh(model,"Hmax",.5,"GeometricOrder","linear");
disp(strcat('Number of Nodes = ','',num2str(length(mesh.Nodes))));
applyBoundaryCondition(model,"neumann", ...
    "Face",1:model.Geometry.NumFaces, ...
    "g",gfun);
% phi = @(xi, eta) [1-xi, 1+xi, 1-eta, 1+eta];
% phi = @(xi, eta) [(1 - xi) * (1 - eta) / 4, (1 + xi) * (1 - eta) / 4,...
%     (1 + xi) * (1 + eta) / 4, (1 - xi) * (1 + eta) / 4];
phi = @(xi, eta) [(1 - xi - eta), xi, eta];
% gaussKernel = fspecial('gaussian', [5 5], 1);
conv_kernel = [];
lambda = 0.1;
tol = 1e-5;
i_sigma = 10;
for i = 1:i_sigma
    eps = 10^(-i);
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
    soln_FE_n_tikh(:,i) = FEMn.B*(FcKinv) + FEMn.ud;
    
    %     soln_FE_s = FEMs.Ks\FEMs.Fs;
    results = solvepde(model,phi);
    
    % extract boundary solution:
    Nf2 = findNodes(mesh,"region","Edge",1:model.Geometry.NumFaces);
    phi_forward_FE_raw(:,i) = results.NodalSolution(Nf2);
    phi_forward_FE_tikh(:,i) = soln_FE_n_tikh(Nf2,i);
    
    soln_FE_raw(:,i) = results.NodalSolution;
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
end

% figure(plot_num)
% % pdemesh(model,"NodeLabels","on")
% % hold on
% plot(mesh.Nodes(1,Nf2),mesh.Nodes(2,Nf2),"ok","MarkerFaceColor","g")
% plot_num = plot_num + 1;


% figure(plot_num)
% pdeplot(model,"XYData",soln_FE_s)
% plot_num = plot_num + 1;

figure(600)
for i = 1:i_sigma
plot(phi_forward_FE_raw(:,i),'DisplayName',strcat('Unregularised, eps =',num2str(10^(-i)))); hold on;
end
xlabel('boundary nodes')
ylabel('\phi(\partial \Omega)')
legend show


figure(601)
for i = 1:i_sigma
plot(phi_forward_FE_tikh(:,i),'DisplayName',strcat('Tikhonov regularised, eps =',num2str(10^(-i)))); hold on;
end
xlabel('boundary nodes')
ylabel('\phi(\partial \Omega)')
legend show

figure(602)
scatter(model.Mesh.Nodes(1,:)',model.Mesh.Nodes(2,:)')

figure(603)
pdeplot(model,"XYData",soln_FE_raw(:,1))
xlabel('x')
ylabel('y')
legend('unregularised')
legend show

figure(604)
pdeplot(model,"XYData",soln_FE_n_tikh(:,1))
xlabel('x')
ylabel('y')
legend('Tikhonov regularised')
legend show
