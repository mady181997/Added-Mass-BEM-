function vertexArea = compute_vertex_areas(F, V)
    % Compute a vertex-based quadrature: each vertex gets one-third the area of each adjacent face.
    numV = size(V,1);
    vertexArea = zeros(numV,1);
    for i = 1:size(F,1)
        idx = F(i, :);
        v1 = V(idx(1), :); v2 = V(idx(2), :); v3 = V(idx(3), :);
        faceArea = 0.5 * norm(cross(v2 - v1, v3 - v1));
        vertexArea(idx) = vertexArea(idx) + faceArea / 3;
    end
end