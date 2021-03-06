function candidate_lines = ...
                findCandidateLines(  strokes_topology,...
                                     intersections,...
                                     pairsInterInter,...
                                     cur_stroke,...
                                     directional_prior,...
                                     cam_param)

global GENERATE_CANDIDATE_PLANES;

%% Initialise
num_intersections = size(cur_stroke.inds_intrsctns_eval,1);

%% Check if there any intersections with previous strokes:

if num_intersections  == 0
    candidate_lines = [];
    return;
end                      

%% First form all the lines with the direction of the direction prior from each of the candidate itersecions:
if ~isempty(directional_prior)    
    
    candidate_lines ...
        = allLinesPassingThroughEachIntersectionVpDirection(...
                intersections,...
                strokes_topology,...
                cur_stroke,...
                directional_prior,...
                cam_param);
else
    candidate_lines = [];
end

disp('allLinesPassingThroughEachIntersectionVpDirection');

%% Enumerate all potential lines formed by pairs of intersections:
   if num_intersections > 1
       candidate_lines = ...
                allLinesFormedByPairsOfIntersections(cur_stroke,...
                                                     strokes_topology,...
                                                     candidate_lines,...
                                                     intersections,...
                                                     cam_param,...
                                                     pairsInterInter);
   end





disp('allLinesFormedByPairsOfIntersections');
if GENERATE_CANDIDATE_PLANES
if (~ismember(cur_stroke.line_group, [1,2,3])) & ~isempty(candidate_lines)
candidate_lines_ = allLinesDominantPlanes(intersections, ...
                                      strokes_topology, ...
                                      cam_param, ...
                                      cur_stroke);
                                      
candidate_lines = [candidate_lines candidate_lines_];
end
end

%% Enumerate all potential lines for lines with prior on the plane:
if (cur_stroke.line_group == 5)
    candidate_lines_ = allLinesPlanePrior(intersections, ...
                                          strokes_topology, ...
                                          cam_param, ...
                                          cur_stroke);

    num_clp = length(candidate_lines_);
    candidate_lines =[candidate_lines candidate_lines_];
end
if isempty(candidate_lines) 
    return;
end


%% Find which nodes can be included in each of the candidate lines:
    
candidate_lines = findAllIntersectionsNearTheLine(cur_stroke, ...
                                            intersections,...
                                            candidate_lines,...
                                            cam_param,...
                                            [strokes_topology(:).depth_assigned],...
                                            strokes_topology);      

disp('findAllIntersectionsNearTheLine');
%% Remove the 3D candidate lines with ubnormally long 3D lengths:
mask = logical(cat(1,strokes_topology(:).depth_assigned));

lns_lngths_3D = cat(1,strokes_topology(mask).length3D);                        

lns_lngths_2D = cat(1,strokes_topology(mask).length2DPrimitive);

if ~isempty(lns_lngths_3D)
    avg_3D_ln_lngth = median(lns_lngths_3D);                        
    ratios3D = cat(1,candidate_lines.length3D)/avg_3D_ln_lngth;
    ratio2D = cur_stroke.length2DPrimitive/median(lns_lngths_2D); 
    lines_keep = ratios3D < 5*ratio2D;
    if sum(lines_keep)
        candidate_lines = candidate_lines(lines_keep);
    end
end 

% fig_num = 10;
% plotCandidateLines(candidate_lines, ...
%                             strokes_topology, ...
%                             cur_stroke,...
%                             intersections, fig_num);
                        
%% Remove the lines with the same set of nodes and the same directions:
    
candidate_lines = ...
         mergeCandidateLines( candidate_lines,...
                              intersections,...
                              cur_stroke,...
                              cam_param,...                              
                              strokes_topology);

disp('done');
end