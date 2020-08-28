function [normals,...
         p_ortho,...
         p_tangent,...
         p_plane] = computePropertiesAtNewIntersection(line_dir,... direction of the stroke where the intersection is added
                                                     candidate_lines, ... %new_stroke
                                                     ind_inter,...
                                                     ind_cnddt_ln, ...
                                                     line_group,...
                                                     intersections, ...
                                                     ind_stroke,...       %new_stroke 
                                                     strokes_topology)
%Compute the normals at the new intersection 'ind_inter'
normals = [];
ind_nrml = 1;

% Max values over all configurations where this intersection is met:
p_ortho = -Inf;
p_tangent = -Inf;
p_plane = -Inf;

inter = [ind_inter; ind_cnddt_ln];

for ind_ln = 1:length(candidate_lines)
    
    for ind_config = 1:length(candidate_lines(ind_ln).configurations)
        configuration = candidate_lines(ind_ln).configurations(ind_config);
        % All the intersections in a candidate line for newly added stroke:
        all_intersections = configuration.inds_intrsctns__mult_cnddts;
        all_intersections_versions = configuration.inds_intrsctns__mult_cnddts_ind;

        mask = ismember(all_intersections, ind_inter);
        if ~sum(mask)
            continue;
        end

        pairs = [all_intersections(mask); all_intersections_versions(mask)];

        if ~(ismember(inter', pairs', 'rows' ))
            continue;
        end

        %Compute the normal between candaite line that has this intersection
        %for new stroke and for old stroke where the intersection is added to.
        normals(ind_nrml,:) = [cross(candidate_lines(ind_ln).dir, line_dir) ind_stroke];
        ind_nrml = ind_nrml+1;

        if line_group == 4
           % line not going towards vanishing point: 
           p_ortho = max(p_ortho, computeScoreOrthogonality(line_dir, candidate_lines(ind_ln).dir));
           p_tangent = max(p_tangent, computeScoreTangential(line_dir, candidate_lines(ind_ln).dir));  
                      
           p_plane = max(p_plane,...
                         checkPlanarityGivenConfiguration(line_dir, ...
                                                          configuration, ...
                                                          intersections, ...
                                                          ind_stroke,...
                                                          strokes_topology,...
                                                          ind_inter));
           %lines_dir -- old stroke, configuration -- new stroke
        end
    end
end
end