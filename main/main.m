clc; clear; close all;
%% Author Madhu Venkata Sri Prudhvi Bikkina

%% Read STL File
[F, V] = readSTL(); % Get faces and vertices from STL file
V = V*0.001; % Convert to meters

%% (Optional) Rotate Vertices
%V = rotate_vertices(V, [1, 0, 0], 90);  % Rotate clockwise (-90 and vice versa) deg about Z-axis
% --- Optional: Check current sphere radius ---
% Assuming the sphere is centered, you can compute its radius as the maximum distance
% from the center of the points:
% center = mean(V, 1);                % Compute the center of the vertices.
% distances = sqrt(sum((V - center).^2, 2));  
% currentRadius = max(distances);
% fprintf('Current sphere radius: %f m\n', currentRadius);
% 
% % --- Determine the scaling factor ---
% desiredRadius = 1.0;  % The radius you want (in meters)
% scaleFactor = desiredRadius / currentRadius;
% fprintf('Scaling factor: %f\n', scaleFactor);
% 
% % --- Scale the vertices ---
% V = V * scaleFactor;

%% Compute Vertex Normals and Vertex Areas
normals = compute_normals(F,V);
vertexArea = compute_vertex_areas(F,V);
fprintf('Total surface area: %f\n', sum(vertexArea));
numV = size(V,1);

%% Precompute Pairwise Differences and Distances
Dx = bsxfun(@minus, V(:,1), V(:,1)'); % Difference in x
Dy = bsxfun(@minus, V(:,2), V(:,2)'); % Difference in y
Dz = bsxfun(@minus, V(:,3), V(:,3)'); % Difference in z
rnorm = sqrt(Dx.^2 + Dy.^2 + Dz.^2); % Pairwise Euclidean distances

%% Assemble the Collocation Matrix (A_coll)
dot_prod = Dx .* repmat(normals(:,1)', numV, 1) + ...
           Dy .* repmat(normals(:,2)', numV, 1) + ...
           Dz .* repmat(normals(:,3)', numV, 1);

A_coll = - ( dot_prod .* repmat(vertexArea', numV, 1) ) ./ (4*pi*(rnorm.^3 + eye(numV)));
A_coll(1:numV+1:end) = 1/2;

%% Precompute Green's Function Matrix
bbox = [max(V(:,1))-min(V(:,1)), max(V(:,2))-min(V(:,2)), max(V(:,3))-min(V(:,3))];
singular_val = 1/(4*pi*mean(bbox));
G = 1./(4*pi*rnorm);
G(1:numV+1:end) = singular_val;

%% Fluid density
rho = 1;
%% Initialize Output Matrices
lambdaMat = zeros(6,6);
phi_all = zeros(numV,6);

%% Loop over the 6 Degrees of Freedom (DOF)
for dof = 1:6
    %% Prescribed Boundary Conditions
    if dof <= 3
        % Translation in x, y, or z direction
        bc = -normals(:, dof);
    else
        % Rotation about axis (dof-3)
        rotAxis = zeros(1,3);
        rotAxis(dof-3) = 1;
        bc = -dot(cross(V, repmat(rotAxis, numV, 1), 2), normals, 2);
    end
    
    %% Solve for potential
    rhs = G * (bc .* vertexArea);
    phi = A_coll \ rhs;
    phi_all(:, dof) = phi;
    
    %% Calculate forces and moments for all directions
    for i = 1:6
        if i <= 3
            % Force calculation
            F_int = -rho * (normals(:,i)' * (phi .* vertexArea));
            lambdaMat(i, dof) = F_int;
        else
            % Moment calculation
            rotAxis = zeros(1,3);
            rotAxis(i-3) = 1;
            r_cross = cross(V, repmat(rotAxis, numV, 1), 2);
            M_int = -rho * sum(dot(r_cross, normals, 2) .* phi .* vertexArea);
            lambdaMat(i, dof) = M_int;
        end
    end
end

%% Ensure matrix symmetry (numerical errors might cause small asymmetries)
lambdaMat = (lambdaMat + lambdaMat')/2;

%% Display the 6x6 Added-Mass Matrix
fprintf('\nAdded-Mass Matrix:\n');
disp(lambdaMat);

%% Visualization: Plot the Surface Potential for each DOF
figure('Position', [100 100 1200 800]);
for dof = 1:6
    subplot(2,3,dof);
    trisurf(F, V(:,1), V(:,2), V(:,3), phi_all(:,dof), 'EdgeColor','none');
    axis equal; colorbar;
    xlabel('x'); ylabel('y'); zlabel('z');
    if dof <= 3
        title(sprintf('Translation DOF %d', dof));
    else
        title(sprintf('Rotation DOF %d', dof-3));
    end
    view(3); lighting phong; camlight headlight;
end

%% --- Helper Function: Rotate Vertices ---
function V_rot = rotate_vertices(V, axis_vec, angle_deg)
    angle_rad = deg2rad(angle_deg);
    axis_vec = axis_vec / norm(axis_vec);
    % Rodrigues' rotation formula
    K = [0, -axis_vec(3), axis_vec(2);
         axis_vec(3), 0, -axis_vec(1);
        -axis_vec(2), axis_vec(1), 0];
    R = eye(3) + sin(angle_rad)*K + (1-cos(angle_rad))*(K*K);
    V_rot = (R * V')';
end