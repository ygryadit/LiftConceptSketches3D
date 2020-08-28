function [strokes_topology,...
          intersections] = assignDepthJointStroke(  intersections,...
                                                    strokes_topology,...
                                                    assign_table,...
                                                    i,...
                                                    cam_param,...
                                                    pairsInterInter)
% assign_table.intrsctns_mlt_hpthss
% assign_table.strks_mlt_hpthss
% assign_table.inds_cnddt_intrsctn           

%%% Intersecting stroke to assign
% ind_intrsctng_strk = assign_table.strks_mlt_hpthss(i);
% ind_intrsctng_strk_pair = intersections(assign_table.intrsctns_mlt_hpthss(i)).strokes_indices == ind_intrsctng_strk;
% 
% % ind_intrsctng_strk = find(intersections(assign_table.intrsctns_mlt_hpthss(i)).strokes_indices == strks_mlt_hpthss(i));
%     
% if (strokes_topology(ind_intrsctng_strk).depth_assigned)
%    return;
% end
% 
% %%% Index of candidate line to select
% 
% inds_cnddt_ln = intersections(assign_table.intrsctns_mlt_hpthss(i)).cnddts3D(assign_table.inds_cnddt_intrsctn(i)).cnddt_lns{ind_intrsctng_strk_pair};
% 
% %%% Indices of valid configurations
% ind_configurations_with_activated_intersection = cell2mat(intersections(assign_table.intrsctns_mlt_hpthss(i)).cnddts3D(assign_table.inds_cnddt_intrsctn(i)).cnfgrtns{ind_intrsctng_strk_pair});
% 
% %%% Select subset of relevant configurations
% candidate_line = strokes_topology(assign_table.strks_mlt_hpthss(i)).candidate_lines(inds_cnddt_ln);
% candidate_line.configurations = candidate_line.configurations(ind_configurations_with_activated_intersection);
% warning('Check that indixing is not messed up');

candidate_line = strokes_topology(assign_table.strks_mlt_hpthss(i)).candidate_lines;

% if length(candidate_line) > 1
   %Find the one with maximum cost.
   costs = zeros(length(candidate_line),1);
   for j = 1:length(candidate_line)
       % inds_non_assigned = find(isempty(cat(1,candidate_line(j).configurations(:).p_full_joint)));
       
       % if ~isempty(inds_non_assigned)
           % for iii = inds_non_assigned(:)
               % try
					% candidate_line(j).configurations(iii).p_full_joint = 0; % the last added stroke has small cost and the branch was not explored
               % catch e
					% rethrow(e);
               % end
           % end
       % end
       max_val = max(cat(1,candidate_line(j).configurations(:).p_full_joint));       
%        if isempty(max_val)
%            max_val =  max(cat(1,candidate_line(j).configurations(:).p_full));
%        end
       costs(j) = max_val;
   end
   [~,ind] = max(costs);
   candidate_line = candidate_line(ind);
% end



%%% Fill in the structure for the current stroke
stroke_joint = strokes_topology(assign_table.strks_mlt_hpthss(i));
stroke_joint.ind = assign_table.strks_mlt_hpthss(i);
UP_TO_LAST = true;
[stroke_joint.inds_intrsctns_eval,...
 stroke_joint.inds_intrsctns_eval_actv,...
 stroke_joint.inds_intrsctns_eval_mltpl_cnddts,...
 stroke_joint.inds_intrsctng_strks_eval,...
 stroke_joint.inds_intrsctng_strks_eval_actv,...
 stroke_joint.inds_intrsctng_strks_eval_mltpl_cnddts] = ...        
    returnIndicesNodesTypes(stroke_joint, ...
                    cat(1, strokes_topology(:).depth_assigned),...
                                        intersections, ...
                                       UP_TO_LAST);


%%% The proper set of configurations will be selcted later:
p_joint_configurations_ = cat(1,candidate_line.configurations(:).p_full_joint);
max_cost_joint = max(cat(1,candidate_line.configurations(:).p_full_joint));


if isempty(candidate_line.configurations)
   disp(''); 
end

try
[strokes_topology,...
 intersections] =...
    assignDepthJointly(strokes_topology, ...
                       intersections, ...
                       cam_param, ...
                       candidate_line, ...
                       stroke_joint, ...
                       p_joint_configurations_,...
                       max_cost_joint,...
                       pairsInterInter); 
catch e
    rethrow(e);
end
% assigned_strokes((end+1):(end+length(vals))) = vals;

%%% At this point depth is assigned to both intersecting strokes
if ~isempty(assign_table.intrsctns_mlt_hpthss(i))
    ind_stroke_1 = intersections(assign_table.intrsctns_mlt_hpthss(i)).strokes_indices(1);
    ind_stroke_2 = intersections(assign_table.intrsctns_mlt_hpthss(i)).strokes_indices(2);
    
    if ~(strokes_topology(ind_stroke_1).depth_assigned) || ...
       ~(strokes_topology(ind_stroke_2).depth_assigned)  
        error('Problem with joint dependent strokes assignment');
    end
    
    intersections(assign_table.intrsctns_mlt_hpthss(i)).cnddts3D = []; %clean 
end

end