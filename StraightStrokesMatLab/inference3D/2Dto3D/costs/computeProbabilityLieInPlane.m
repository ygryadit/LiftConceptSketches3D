function p = computeProbabilityLieInPlane(normals, candidate_line_dir)
    global sigma_costs;
    
    if isempty(normals)
        p = 0;
        return;
    end
    
    candidate_line_dir = candidate_line_dir/norm(candidate_line_dir);
    candidate_line_dir = repmat(candidate_line_dir, size(normals,1), 1);

    cos_angles = abs(dot(normals, candidate_line_dir, 2));
    cos_angle = min(cos_angles);
%     sigma = cos(pi/2 - pi/24);

    p = exp(-(cos_angle).^2/(2*sigma_costs^2));

end