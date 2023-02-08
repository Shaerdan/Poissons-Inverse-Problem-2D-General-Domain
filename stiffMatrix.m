function [K] = stiffMatrix(dt)
% Define number of nodes and number of elements
nnodes = length(dt.Points);
nelements = length(dt.ConnectivityList);

% Initialize global stiffness matrix K
K = sparse(nnodes, nnodes);

% Loop over each triangular element
for e = 1:nelements
    % Get indices of nodes of the current element
    nodes = dt.ConnectivityList(e,:);
    
    
    % Get the x, y coordinates of the nodes of the current element
    x1 = dt.Points(nodes(1),1);
    y1 = dt.Points(nodes(1),2);
    x2 = dt.Points(nodes(2),1);
    y2 = dt.Points(nodes(2),2);
    x3 = dt.Points(nodes(3),1);
    y3 = dt.Points(nodes(3),2);
    
    % Compute the area of the current element
    Area = abs(0.5 * (x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)));
    % Define natural coordinates xi and eta
    xi = [1, x1, x2, x3];
    eta = [1, y1, y2, y3];
    % Define finite element basis functions for linear elements
    N = [1-xi-eta, xi, eta];
    dNdxi = [-1, -1; 1, 0; 0, 1];
    % Evaluate gradient of basis functions at nodes
    J = [[x1, x2, x3]*dNdxi(:,1); [y1, y2, y3]*dNdxi(:,1)];
    dNdx = dNdxi * inv(J);
    % Compute element stiffness matrix
    Ke = zeros(3, 3);
    for i = 1:3
        for j = 1:3
            Ke(i,j) = sum(dNdx(i,:).*dNdx(j,:),2) * Area;
        end
    end
    
    % Assemble element stiffness matrix into global stiffness matrix
    K(nodes, nodes) = K(nodes, nodes) + Ke;
    
end
end