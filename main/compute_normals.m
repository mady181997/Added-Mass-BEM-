function normals = compute_normals(F, V)
    % Compute vertex normals by averaging the normals of adjacent faces.
    numV = size(V,1);
    normals = zeros(numV, 3);
    for i = 1:size(F,1)
        idx = F(i, :);
        v1 = V(idx(1), :); v2 = V(idx(2), :); v3 = V(idx(3), :);
        fn = cross(v2 - v1, v3 - v1);
        fn = fn / norm(fn);
        normals(idx, :) = normals(idx, :) + repmat(fn, 3, 1);
    end
    normals = normals ./ vecnorm(normals, 2, 2);
end