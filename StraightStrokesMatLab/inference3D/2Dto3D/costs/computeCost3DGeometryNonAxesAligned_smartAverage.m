function p_directional = computeCost3DGeometryNonAxesAligned_smartAverage(...
                                                candidate_line,...
                                                pairsInterInter,...
                                                strokes_topology,...
                                                intersections,...
                                                cur_stroke)
global DIRECTIONAL_NON_AXIS;
global max_cost_threshold;
num_configs = length(candidate_line.configurations);
p_directional = zeros(num_configs,1);                                            
                                            
if ~DIRECTIONAL_NON_AXIS
    return;
end

% First all the intersections in all configurations with assigned lines:
mask = arrayfun(@(x) ~isempty(x.inds_intrsctns__assigned), candidate_line.configurations);
inds_intrsctns__assigned_all = ...
    unique(cat(2,candidate_line.configurations(mask).inds_intrsctns__assigned));

strokes_at_inter = cat(1,intersections(inds_intrsctns__assigned_all).strokes_indices);
strokes_at_inter = strokes_at_inter(strokes_at_inter ~= cur_stroke.ind);
        

% Second all the intersections in all configurations with lines with
% multipe caniddates:
mask_mult_cnddts =  arrayfun(@(x) ~isempty(x.inds_intrsctns__mult_cnddts), candidate_line.configurations);

inds_intrsctns__mult_cnddts = ...
    [cat(2,candidate_line.configurations(mask_mult_cnddts ).inds_intrsctns__mult_cnddts);...
     cat(2,candidate_line.configurations(mask_mult_cnddts ).inds_intrsctns__mult_cnddts_ind)]';
 
inds_intrsctns__mult_cnddts = unique(inds_intrsctns__mult_cnddts, 'rows');


%% Precompute ortho and tangent for each intersection:
%  I. Go only over the intersectiosn that are actually in the line:                                             
num_intrsctns_1 = length(inds_intrsctns__assigned_all);
num_intrsctns_2 = size(inds_intrsctns__mult_cnddts,1);

p_ortho = zeros(num_intrsctns_1+num_intrsctns_2,1);
p_tangent = zeros(num_intrsctns_1+num_intrsctns_2,1);

line_dir_1 = candidate_line.dir;

for i = 1:num_intrsctns_1
    
    line_dir_2 = strokes_topology(strokes_at_inter(i)).direction_vec;
    p_ortho(i) = computeScoreOrthogonality(line_dir_1, line_dir_2);
    p_tangent(i) = computeScoreTangential(line_dir_1, line_dir_2);   
end

% II. Go over the intersections with strokes with multiple candidates:
if ~isempty(inds_intrsctns__mult_cnddts)
    strokes_at_inter2 = cat(1,intersections(inds_intrsctns__mult_cnddts(:,1)).strokes_indices);
    inds_swap = find(strokes_at_inter2(:,1)) == cur_stroke.ind;
    strokes_at_inter2(inds_swap,1) = strokes_at_inter2(inds_swap,2);
    strokes_at_inter2(inds_swap,2) = cur_stroke.ind;
    
    strokes_at_inter2 = strokes_at_inter2(:,1);

    for i = num_intrsctns_1+(1:num_intrsctns_2)
        %Preinitialise
        p_ortho(i) = -Inf;
        p_tangent(i) = -Inf;

        %Find candidate lines indices.
        i2 = i - num_intrsctns_1;
        mask_pair = intersections(inds_intrsctns__mult_cnddts(i2,1)).strokes_indices == strokes_at_inter2(i2);
        try
            inds_cnddt_lns = intersections(inds_intrsctns__mult_cnddts(i2,1)).cnddts3D(inds_intrsctns__mult_cnddts(i2,2)).cnddt_lns{mask_pair};
        catch e
            rethrow(e);
        end
        
        for icl = inds_cnddt_lns
            try
                line_dir_2 = strokes_topology(strokes_at_inter2(i2)).candidate_lines(icl).dir;
            catch e
                rethrow(e);
            end
            p_ortho(i) = max(p_ortho(i), computeScoreOrthogonality(line_dir_1, line_dir_2));
            p_tangent(i) = max(p_tangent(i), computeScoreTangential(line_dir_1, line_dir_2));   
        end
    end
end

%% Planarity for each configration
P_plane = checkLinePlanarity(candidate_line, ...
                                      intersections, ...
                                      strokes_topology, ...
                                      cur_stroke.ind);

%% Go over configurations and select maximum:
num_good_intersections = zeros(num_configs,1);
p_geom = zeros(num_configs,1);

for i = 1:num_configs
    configuration = candidate_line.configurations(i);
    p_directional(i) = -Inf;
    
    ind1 = find(ismember(inds_intrsctns__assigned_all, configuration.inds_intrsctns));
    
    cngrtn_inds_intrsctns__mult_cnddts = [candidate_line.configurations(i).inds_intrsctns__mult_cnddts;...
                                   candidate_line.configurations(i).inds_intrsctns__mult_cnddts_ind]';
    if ~isempty(inds_intrsctns__mult_cnddts) & ~isempty(cngrtn_inds_intrsctns__mult_cnddts)
        ind2 = find(ismember(inds_intrsctns__mult_cnddts, cngrtn_inds_intrsctns__mult_cnddts, 'rows'));
    else
        ind2 = [];
    end
    clear('p_vals');
    p_vals(1,:) = [p_ortho(ind1)'  p_ortho(ind2+ num_intrsctns_1)'];
    p_vals(2,:) = [p_tangent(ind1)' p_tangent(ind2+ num_intrsctns_1)'];
 
    p_vals = max(p_vals,[],2);
    p_vals = p_vals(p_vals > max_cost_threshold);
    num_good_intersections(i) = length(p_vals);
    
    p_geom(i) = sum(p_vals);
    
%     for ind_intrsctn = configuration.inds_intrsctns
%         
%         % Find_intersection
%         ind = find(ismember(inds_intrsctns__assigned_all, ind_intrsctn)); 
%         
%         if ~isempty(ind)
%             p_ortho_ = p_ortho(ind);
%             p_tangent_ = p_tangent(ind);
%         else
%             ind_ = configuration.inds_intrsctns__mult_cnddts_ind(configuration.inds_intrsctns__mult_cnddts == ind_intrsctn);
%             
%             ind = find(ismember(inds_intrsctns__mult_cnddts, [ind_intrsctn, ind_], 'rows'));
%             p_ortho_ = p_ortho(ind+num_intrsctns_1);
%             p_tangent_ = p_tangent(ind+num_intrsctns_1);
%         end
%         
%         p_directional(i) = max([p_directional(i),P_plane(i),p_ortho_,p_tangent_]); 
%         
%     end
end

if max(num_good_intersections)~=0
p_geom = p_geom./max(num_good_intersections);
p_directional = 0.5*(P_plane + p_geom);
else
    p_directional = P_plane;
end





% %% Go over configurations and find groups of intersections:
% 
% for i = 1:num_configs
%     configuration = candidate_line.configurations(i);
%     p_directional(i) = -Inf;
%     
%     for ind_intrsctn = configuration.inds_intrsctns
%         
%         % Planarity
%         inds_cls_intrsctns = findAllIntersectionsInNeighborhood(pairsInterInter,...
%                                         ind_intrsctn,...
%                                         configuration.inds_intrsctns);
%         
%         try
%         dirs = findDirectionsInIntersectionsGroup(inds_cls_intrsctns,...
%                                             configuration.inds_intrsctns__assigned,...
%                                             configuration.inds_intrsctns__mult_cnddts,...
%                                             configuration.inds_intrsctns__mult_cnddts_ind,...
%                                             intersections,...
%                                             strokes_topology,...
%                                             cur_stroke.ind);
%         catch e
%             rethrow(e);
%         end
%         
%         if size(dirs,1) >1
%             normals = computeNormalsIntersection(dirs);
%             
%             p_plane = ...
%                  computeProbabilityLieInPlane(normals, line_dir_1);
%         else
%             p_plane = 0.0;
%         end
%         
%         % Find_intersection
%         ind = find(ismember(inds_intrsctns__assigned_all, ind_intrsctn)); 
%         
%         if ~isempty(ind)
%             p_ortho_ = p_ortho(ind);
%             p_tangent_ = p_tangent(ind);
%             
%         else
%             ind_ = configuration.inds_intrsctns__mult_cnddts_ind(configuration.inds_intrsctns__mult_cnddts == ind_intrsctn);
%             
%             ind = find(ismember(inds_intrsctns__mult_cnddts, [ind_intrsctn, ind_], 'rows'));
%             p_ortho_ = p_ortho(ind+num_intrsctns_1);
%             p_tangent_ = p_tangent(ind+num_intrsctns_1);
%         end
%         
%         p_directional(i) = max([p_directional(i),p_plane,p_ortho_,p_tangent_]);        
% %         fprintf('ind_intrsctn %d, p_plane %.2f, p_ortho_ %.2f, p_tangent_ %.2f, p_directional(%d) %.2f\n', ...
% %             ind_intrsctn, p_plane,p_ortho_,p_tangent_, i,  p_directional(i));
%         
%     end
%     
%     
% end



return;

%     num_intersections = length(inlines_intesrections_indices);
% 
%     strokes_mult_hypothesis = cur_stroke.inds_intrsctng_strks_eval_mltpl_cnddts;
%     
%     
%     cur_stroke_ind = cur_stroke.ind;
%     
%     p_intersection = zeros(num_intersections,1);
%     
%     for i = 1:num_intersections
%       
%         int_ind = inlines_intesrections_indices(i);
%         
%         % Find a local group of intersections to find all the strokes
%         % meeting at that point:
%         mask_inds_cls_intrsctns = (pairsInterInter(:,1) == int_ind);
%         inds_cls_intrsctns = pairsInterInter(mask_inds_cls_intrsctns,2);
%         mask_inds_cls_intrsctns = (pairsInterInter(:,2) == int_ind);
%         inds_cls_intrsctns = [inds_cls_intrsctns; pairsInterInter(mask_inds_cls_intrsctns,1)];
%         
%         % Keep here only intersections with assigned strokes (no strokes with multiple hypothesis):
%         inds_cls_intrsctns = inds_cls_intrsctns(cat(1,intersections(inds_cls_intrsctns).is_active)==1);
%         
% %         inds_cls_intrsctns = [inds_cls_intrsctns; int_ind];
%         
%         
%         % Find all the assigned strokes at the group:
%         strokes_at_inter = cat(1,intersections(inds_cls_intrsctns).strokes_indices);
%         strokes_at_inter = setdiff(strokes_at_inter(:), cur_stroke_ind, 'stable');
%           
%         [mask_inter_mult_hypoth, inds_strokes] = ismember(strokes_at_inter, strokes_mult_hypothesis);
%         intersecting_strokes_dirvecs = cat(1,strokes_topology(strokes_at_inter(~mask_inter_mult_hypoth)).direction_vec);
% %         intersecting_strokes_dirvecs = cat(1,strokes_topology(strokes_at_inter).direction_vec);
%         
%         strokes_at_intersection = cat(1,intersections(int_ind).strokes_indices);
%         
%         stroke_at_cur_inter = setdiff(strokes_at_intersection(:), cur_stroke_ind, 'stable');
%         
%         if strokes_topology(stroke_at_cur_inter).depth_assigned
%             intersecting_strokes_dirvecs(end+1,:) = strokes_topology(stroke_at_cur_inter).direction_vec;
%         else
%             
%             ind_intrsctn_cnddt_crdnts = inds_intrsctns__mult_cnddts_ind(inds_intrsctns__mult_cnddts == int_ind);
%             mhs = find(intersections(int_ind).strokes_indices == stroke_at_cur_inter);
%             ln_hypoth_num = intersections(int_ind).cnddts3D(ind_intrsctn_cnddt_crdnts).cnddt_lns{mhs};
% %             hypoth_num = intersections%ind_hypothesis_in_stroke(inter_mult_hypothesis == int_ind);
%             
%             for lhmj = ln_hypoth_num
%                 intersecting_strokes_dirvecs(end+1,:) = strokes_topology(stroke_at_cur_inter).candidate_lines(lhmj).dir;
%             end
%         end
%         
%         
%         %Find the direction of the stroke accordig to consistent hypothesis: 
% %         if sum(mask_inter_mult_hypoth)
% %             curr_str_strokes_mult_hypoths = strokes_at_inter(mask_inter_mult_hypoth);
% %             inds_strokes = inds_strokes(mask_inter_mult_hypoth);
% %             
% %             for jj = 1:length(curr_str_strokes_mult_hypoths)
% %                 hsjj = ind_hypothesis_in_stroke(inds_strokes(jj))
% %                 sjj = curr_str_strokes_mult_hypoths(jj)
% %                 if ~isempty(strokes_topology(sjj).candidate_lines)
% %                     intersecting_strokes_dirvecs(end+1,:) = strokes_topology(sjj).candidate_lines(hsjj).dir;
% %                 end
% %             end
% %         end
% %         figure(49);   
% %         hold off;
% %         plot3([0,line_dir(1)], [0,line_dir(2)],[0,line_dir(3)],'r'); hold on;
% %         
%         % Orthogonality:
%         for j = 1:size(intersecting_strokes_dirvecs,1)
%             
% %              plot3( [0,intersecting_strokes_dirvecs(j,1)],...
% %                     [0,intersecting_strokes_dirvecs(j,2)],...
% %                     [0,intersecting_strokes_dirvecs(j,3)],'g');
% %             
%             cos_angle = abs(dot(line_dir, intersecting_strokes_dirvecs(j,:)));
%             sigma = cos(pi/2 - pi/36);
% 
%             p_intersection(i) = max(exp(-(cos_angle).^2/(2*sigma^2)),p_intersection(i));
% %             p = max(p, p_intersection);
% %             p = probabilityAorB(p, p_intersection);
%         end
%         
%         % Belongs to a plane
%         if size(intersecting_strokes_dirvecs,1) >1
%             normals = computeNormalsIntersection(intersecting_strokes_dirvecs);
%             
%             p_plane = ...
%                  computeProbabilityLieInPlane(normals, line_dir);
%     %         p = probabilityAorB(p, p_plane);    
%             p_intersection(i) = probabilityAorB(p_intersection(i), p_plane);    
%         end
%     end
%     
%     
% %     for i = reshape(ind_snap_intersection_active_i, 1, [])
% %         ind_intersecting_strokes = nodes_grouped_active_incident_strokes(i,:);
% %         direvecs_i = cat(1,wires_strokes(ind_intersecting_strokes).direction_vec);
% %         
% %         normals = computeNormalsIntersection(direvecs_i);
% %         
% %         p_plane = ...
% %              computeProbabilityLieInPlane(normals, line_dir);
% %         p = probabilityAorB(p, p_plane);    
% %     end      
%         
% %         p = mean(p_intersection);
%         p = max(p_intersection);
end






function inds_cls_intrsctns = findAllIntersectionsInNeighborhood(pairsInterInter, int_ind, inds_intrsctns_ln)

% Find a local group of intersections to find all the strokes
% meeting at that point:
mask_inds_cls_intrsctns = (pairsInterInter(:,1) == int_ind);
inds_cls_intrsctns = pairsInterInter(mask_inds_cls_intrsctns,2);
mask_inds_cls_intrsctns = (pairsInterInter(:,2) == int_ind);
inds_cls_intrsctns = [inds_cls_intrsctns; pairsInterInter(mask_inds_cls_intrsctns,1)];
inds_cls_intrsctns = intersect(inds_cls_intrsctns, inds_intrsctns_ln);
inds_cls_intrsctns = [inds_cls_intrsctns', int_ind];
end






% function p = computeProbabilityLieInPlane(normals, candidate_line_dir)
%     if isempty(normals)
%         p = 0;
%         return;
%     end
%     candidate_line_dir = candidate_line_dir/norm(candidate_line_dir);
% 
%     candidate_line_dir = repmat(candidate_line_dir, size(normals,1), 1);
% 
%     cos_angles = abs(dot(normals, candidate_line_dir, 2));
%     cos_angle = min(cos_angles);
%     sigma = cos(pi/2 - pi/12);
% 
%     p = exp(-(cos_angle).^2/(2*sigma^2));
% 
% end