function [strokes_topology, intersections] = ...
                checkDepthAssignemnt(strokes_topology, ...
                                     intersections,...
                                     cur_stroke,...
                                     candidate_lines,...
                                     pairsInterInter,...
                                     cam_param,...                                     
                                     do_assign_depth)


costs = zeros(length(candidate_lines),1);
if isempty(candidate_lines)
    % The case when for some stroke there is no plausible explantion
    % left.
    strokes_topology(cur_stroke.ind).depth_assigned = true;
    strokes_topology(cur_stroke.ind).candidate_lines = [];
    strokes_topology(cur_stroke.ind).num_candidate_lines = 0;
    return;
end


for ii = 1:length(candidate_lines)
    % for i = 1:length(candidate_lines(ii).configurations)
       % if isempty(candidate_lines(ii).configurations(i).p_full_joint) 
           % candidate_lines(ii).configurations(i).p_full_joint = ...
               % candidate_lines(ii).configurations(i).p_full;
       % end
    % end
    costs(ii) = max(cat(1, candidate_lines(ii).configurations(:).p_full_joint));
end

%     costs = cat(1,       candidate_lines(:).max_cost);               
[max_cost,ind_line_most_probable] = max(costs);
candidate_line = candidate_lines(ind_line_most_probable);
        
global DELAY_ASSIGN
    num_candidate_lines = length(candidate_lines);
% if DELAY_ASSIGN 
% 
%     if doNotAssignDepthValueToLine(cur_stroke.score,...
%           strokes_topology(cur_stroke.ind).confidence,...
%           num_candidate_lines,...
%           candidate_lines) || ~do_assign_depth
% 
%           fprintf('Not assigned %d: max_cost = %.3f, confidence = %.3f\n',...
%               cur_stroke.ind, cur_stroke.score, strokes_topology(cur_stroke.ind).confidence)
%           return;
%     end
% 
% end
%     
    %% Assign depth value and remove axilary edges from the graph.
    


    try
        p_confgrtns_ln_mst_prbbl = cat(1, candidate_line.configurations(:).p_full_joint);
    catch
        disp('');
    end
    disp(cur_stroke.score);
    disp(max_cost);
    [strokes_topology, intersections] = assignDepthJointly(strokes_topology,...
                                                               intersections,...
                                                               cam_param,...
                                                               candidate_line,...
                                                               cur_stroke,...
                                                               p_confgrtns_ln_mst_prbbl,...
                                                               max_cost,...
                                                               pairsInterInter); 

end