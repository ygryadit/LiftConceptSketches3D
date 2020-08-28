function p_tangent = computeScoreTangential(line_dir_1, line_dir_2)
    global sigma_aligned;
    cos_angle = abs(dot(line_dir_1, line_dir_2));
%     sigma = cos(pi/2 - pi/24);
    sigma = sind(0.75);
    %p_tangent = exp(-(1.0 - cos_angle.^2)/(2.0*sigma_aligned^2));
    p_tangent = exp(-((1.0 - cos_angle).^2)/(2.0*sigma_aligned^2));
end