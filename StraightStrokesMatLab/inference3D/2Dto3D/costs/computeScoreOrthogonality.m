function p_ortho = computeScoreOrthogonality(lines_dir_1, line_dir_2)
    try
        lines_dir_1 = normilizeDirections(lines_dir_1);
    catch e
        rethrow(e);
    end
    
    num_lines = size(lines_dir_1,1);
    line_dir_2 = repmat(line_dir_2, num_lines, 1);
    
    
    global sigma_costs;
    
    cos_angle = abs(dot(lines_dir_1, line_dir_2, 2));
    %sigma = cos(pi/2 - pi/36);
%     sigma = cos(pi/2 - pi/24);
    
    p_ortho = exp(-(cos_angle).^2/(2*sigma_costs^2));
end

function lines_dirs = normilizeDirections(lines_dirs)
    lengths     = repmat(sqrt(sum(lines_dirs.^2, 2)), 1, 3);
    lines_dirs  = lines_dirs./lengths ;
end