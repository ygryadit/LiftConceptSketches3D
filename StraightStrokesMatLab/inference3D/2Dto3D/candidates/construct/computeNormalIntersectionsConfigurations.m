function candidate_line  = computeNormalIntersectionsConfigurations(candidate_line,...
                                                strokes_topology,...
                                                intersections,...
                                                cur_stroke)
num_configs = length(candidate_line.configurations);
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

%% Precompute normals at each intersection:
num_intrsctns_1 = length(inds_intrsctns__assigned_all);
num_intrsctns_2 = size(inds_intrsctns__mult_cnddts,1);

normals = cell(num_intrsctns_1+num_intrsctns_2,1);


line_dir_1 = candidate_line.dir;
%  I. Go only over the intersectiosn that are actually in the line:   
for i = 1:num_intrsctns_1
     line_dir_2 = strokes_topology(strokes_at_inter(i)).direction_vec;
     normals{i} = [cross(line_dir_1, line_dir_2) strokes_at_inter(i)];
end

j = i;

% II. Go over the intersections with strokes with multiple candidates:
if ~isempty(inds_intrsctns__mult_cnddts)
    strokes_at_inter2 = cat(1,intersections(inds_intrsctns__mult_cnddts(:,1)).strokes_indices);
    inds_swap = find(strokes_at_inter2(:,1)) == cur_stroke.ind;
    strokes_at_inter2(inds_swap,1) = strokes_at_inter2(inds_swap,2);
    strokes_at_inter2(inds_swap,2) = cur_stroke.ind;
    
    strokes_at_inter2 = strokes_at_inter2(:,1);

    for i = num_intrsctns_1+(1:num_intrsctns_2)
           i2 = i - num_intrsctns_1;
        mask_pair = intersections(inds_intrsctns__mult_cnddts(i2,1)).strokes_indices == strokes_at_inter2(i2);
        try
            inds_cnddt_lns = intersections(inds_intrsctns__mult_cnddts(i2,1)).cnddts3D(inds_intrsctns__mult_cnddts(i2,2)).cnddt_lns{mask_pair};
        catch e
            rethrow(e);
        end
        j = 1;
        for icl = inds_cnddt_lns
            try
                line_dir_2 = strokes_topology(strokes_at_inter2(i2)).candidate_lines(icl).dir;
            catch e
                rethrow(e);
            end
            j = j+1;
            normals{i}(j,:) = [cross(line_dir_1, line_dir_2) strokes_at_inter2(i2)];
        end
    end
end

%% Go over the configurations and assign normals
for i = 1:num_configs
     configuration = candidate_line.configurations(i);
     for ind_intrsctn = configuration.inds_intrsctns
         ind = find(ismember(inds_intrsctns__assigned_all, ind_intrsctn)); 
         if ~isempty(ind)
             if norm(normals{ind}(1,1:3)) > 0.5
                configuration.planes_normals(end+1,:) = normals{ind};
             end
         else
             
             ii = ismember(configuration.inds_intrsctns__mult_cnddts, ind_intrsctn);
             
             ind = find(ismember(inds_intrsctns__mult_cnddts,...
                        [ind_intrsctn configuration.inds_intrsctns__mult_cnddts_ind(ii)], 'rows')); 
                    
             for j =1:size(normals{ind+num_intrsctns_1},1)
                 if norm(normals{ind+num_intrsctns_1}(j,1:3)) > 0.5
                    configuration.planes_normals(end+1,:) = normals{ind+num_intrsctns_1}(j,:);
                 end
             end
         end
     end
     
     candidate_line.configurations(i) = configuration;
end
end