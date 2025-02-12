function [faces, vertices] = readSTL()
    % This function reads an STL file selected by the user and returns the faces and vertices.
    % Output: 
    %   faces: The connectivity list (triangle faces).
    %   vertices: The coordinates of the vertices.

    % Open a file selection dialog to choose the STL file
    [fileName, filePath] = uigetfile('*.stl', 'Select an STL file');

    % Check if the user canceled the selection
    if isequal(fileName, 0)
        disp('User canceled the file selection.');
        return;
    end

    % Construct the full path of the selected file
    stlFile = fullfile(filePath, fileName);

    % Read the STL file (returns a structure in newer MATLAB versions)
    stlData = stlread(stlFile);

    % Extract faces and vertices
    faces = stlData.ConnectivityList; % Triangle faces
    vertices = stlData.Points; % Vertex coordinates
   %%  vertices = vertices * 0.001;  % Convert to meters (example scaling)

    % Plot the extracted boundary surface
figure;
trisurf(faces, vertices(:,1), vertices(:,2), vertices(:,3), ...
    'FaceColor', 'cyan', 'EdgeColor', 'k');
title('Extracted Boundary Surface');
axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
end
