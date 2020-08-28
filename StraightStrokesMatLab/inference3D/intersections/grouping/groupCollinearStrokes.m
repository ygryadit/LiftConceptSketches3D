% function groupCollinearStrokes()
% 
% Groups collinear stokes, where they overlap for at least 85 percent of
% length of the shortest stroke.
% Goes over the strokes in the order they are drawn and checks if the
% grouping can be performed.

function [intersections, strokes_topology] = ...
            groupCollinearStrokes(lines_pairs,...
                               intersections,...
                               strokes_topology,...
                               img)
                           
    debug_merge = false;                           
                            
    lines_pairs = rearrangePairs(lines_pairs);
    
    strokes_to_remove = [];
    
    
    if debug_merge
        figure(12);        
        hold off;
        imshow(img);
        hold on;
    end
    
    for i = 1:length(strokes_topology)
        if i == 128
           disp(''); 
        end
        if strokes_topology(i).primitive_type ~= 0
           continue; 
        end
        
        lines_pairs = rearrangePairs(lines_pairs);
    
          
        inds = find((lines_pairs(:,1) == i));
        ind_ln_merge = lines_pairs(inds,2);
        ind_ln_merge = setdiff(ind_ln_merge, i);
        ind_ln_merge = sort(ind_ln_merge, 'ascend');
        
        
        if ~isempty(ind_ln_merge) 
            
            if debug_merge
                disp('-------------');
            end
            ind_strk_1 = i;
            while ~isempty(ind_ln_merge)
%             for kk = 1:length(ind_ln_merge)%length(ind_ln_merge):-1:1
                
                ind_strk_2 = ind_ln_merge(1);
%                 fprintf('ind_strk_1 = %d, ind_strk_2 = %d\n', ind_strk_1, ind_strk_2);
       
%                 if (ind_strk_1 == 18) & (ind_strk_2 == 21)
%                     disp('');
%                 end
%                 if (ind_strk_1 == 18) & (ind_strk_2 == 31)
%                     disp('');
%                 end
%                     
                vec1 = strokes_topology(ind_strk_1).primitive_geom;
                vec2 = strokes_topology(ind_strk_2).primitive_geom;

                [begin, end_coordinate, vec1, vec2, dist] = findEndpointsCommonRegion(vec1, vec2);
                
%                 [length_commom, length_s1, length_s2, s_new_points2D, dist] = ...
%                     findEndpointsCommonRegionPolyPolyLines(strokes_topology(ind_strk_1).points2D, ...
%                                                       strokes_topology(ind_strk_2).points2D,...
%                                                       strokes_topology(ind_strk_2).accuracy_radius);
                                                  

                if debug_merge
                       figure(12);        
                    hold off;
                    imshow(img);
                    hold on;
                    plot([strokes_topology(ind_strk_1).points2D(:).x],...
                        [strokes_topology(ind_strk_1).points2D(:).y], 'b');        
                    plot([strokes_topology(ind_strk_2).points2D(:).x], ...
                        [strokes_topology(ind_strk_2).points2D(:).y], 'g');
                    
                   
                end



                mergeSameGroups = ...
                    strokes_topology(ind_strk_1).line_group == ...
                    strokes_topology(ind_strk_2).line_group;
                
                
                %% Check thac collinear:
                dir1 = vec1([2,4]) - vec1([1,3]);
                dir2 = vec2([2,4]) - vec2([1,3]);
                dir1 = dir1./norm(dir1);
                dir2 = dir2./norm(dir2);
                cos_dirs = dot(dir1, dir2, 2);
                are_collinear = abs(cos_dirs) > cos(pi/36) & mergeSameGroups;

                if ~are_collinear 
                    ind_inter = find( (intersections.strokes_indices(:,1) == ind_strk_2 & intersections.strokes_indices(:,2) == ind_strk_1) | ...
                      (intersections.strokes_indices(:,2) == ind_strk_2 & intersections.strokes_indices(:,1) == ind_strk_1) );
                    intersections.collinear(ind_inter) = false;
%                         fprintf('Changed to not collinear: ');
%                         disp(ind_inter);
                end
                
                %%
                areContinious_ = areContinious(strokes_topology, ...
                                        ind_strk_1,...
                                        ind_strk_2);
                                    
%                 indcs_intrsctng_strks = strokes_topology(ind_strk_1).indcs_intrsctng_strks(strokes_topology(ind_strk_1).indcs_intrsctng_strks > ind_strk_1);
%                 indcs_intrsctng_strks = sort(indcs_intrsctng_strks, 'ascend');
%                 if ~isempty(indcs_intrsctng_strks)
%                     n = min(2, length(indcs_intrsctng_strks));
%                     areContinious_ = ismember(ind_strk_2, indcs_intrsctng_strks(1:n));
%                 else
%                     areContinious_  =false;
%                 end
                
        
                if mergeSameGroups &  are_collinear &...%checkMerge(length_commom, length_s1, length_s2,strokes_topology,ind_strk_1, ind_strk_2, dist)
                    checkMerge(begin, end_coordinate,strokes_topology,ind_strk_1, ind_strk_2, dist, areContinious_)
                    
                    % Merge
                    line = merge(vec1, vec2, ...
                                strokes_topology(ind_strk_1).mean_pressure,...
                                strokes_topology(ind_strk_2).mean_pressure);

%                     axis equal
%                     plot(line([1,2]), line([3,4]), 'Color', 'r', 'LineWidth',2);

                    % Write into the first stroke, 
                    strokes_topology(ind_strk_1).primitive_geom = line;

%                     strokes_topology(ind_strk_1).points2D = ...
%                         [strokes_topology(ind_strk_1).points2D; ...
%                          strokes_topology(ind_strk_2).points2D];
                    strokes_topology(ind_strk_1).points2DOriginal = ...
                        strokes_topology(ind_strk_1).points2D;
                    
                    strokes_topology(ind_strk_1).points2D = ...
                        approximateMergedStrokesSingleLine( ...
                            strokes_topology(ind_strk_1).points2D,...
                            strokes_topology(ind_strk_2).points2D);
                    
                    
                        
                        
                    if debug_merge
                        plot([strokes_topology(ind_strk_1).points2D(:).x],...
                            [strokes_topology(ind_strk_1).points2D(:).y], 'c');        
                        fprintf("Merged %d and %d\n", ind_strk_1, ind_strk_2);
                    end
                
                    strokes_topology(ind_strk_1).mean_pressure = ...
                        max([strokes_topology(ind_strk_1).mean_pressure,...
                             strokes_topology(ind_strk_2).mean_pressure]);

                    strokes_topology(ind_strk_1).length2DFull = norm(line([2,4])- line([1,3]));
                    strokes_topology(ind_strk_1).length2DPrimitive = norm(line([2,4])- line([1,3]));
                    strokes_topology(ind_strk_1).accuracy_radius = max( strokes_topology(ind_strk_1).accuracy_radius,...
                                                                        strokes_topology(ind_strk_2).accuracy_radius) + 0.5*dist;

                    merged_with =    sort(unique([strokes_topology(ind_strk_1).merged_with ...
                                                                       strokes_topology(ind_strk_2).merged_with ...
                                                                       ind_strk_1...
                                                                       ind_strk_2]));
                                                                   
                    strokes_topology(ind_strk_1).merged_with = merged_with;
                    strokes_topology(ind_strk_2).merged_with = merged_with;
                                                                   
                    % Change the strokes indices in the intersections
                    if ~isfield(intersections, 'strokes_indices_original') || isempty(intersections.strokes_indices_original)
                        intersections.strokes_indices_original = intersections.strokes_indices;
                    end
                    
                    ind1 = find(intersections.strokes_indices(:,1) == ind_strk_2 & intersections.strokes_indices(:,2) ~= ind_strk_1);
                    intersections.strokes_indices(ind1,1) = ind_strk_1;
                    ind2 = find(intersections.strokes_indices(:,2) == ind_strk_2 & intersections.strokes_indices(:,1) ~= ind_strk_1);
                    intersections.strokes_indices(ind2,2) = ind_strk_1;
                    
                    
                    % Remove the stroke
                    strokes_to_remove(end+1) = ind_strk_2;
                    % Reindex the strokes indices in lines_pairs
                    ind1 = find(lines_pairs(:,1) == ind_strk_2);
                    lines_pairs(ind1,1) = ind_strk_1;
                    
                    ind1 = find(lines_pairs(:,2) == ind_strk_2);
                    lines_pairs(ind1,1) = ind_strk_1;
                    
                    lines_pairs = rearrangePairs(lines_pairs);
                    lines_pairs = setdiff(lines_pairs, [ind_strk_1 ind_strk_2], 'rows');
                    
                    inds = find((lines_pairs(:,1) == i));
                    ind_ln_merge = lines_pairs(inds,2);
%                     ind_ln_merge = [];
                    
                
                else
                    ind_ln_merge = ind_ln_merge(2:end);
                    lines_pairs = setdiff(lines_pairs, [ind_strk_1 ind_strk_2], 'rows');
%                     disp('not grouped')   
                end
                
                
                
            end
        end
    end
    for i = 1:length(strokes_to_remove)
        strokes_topology(strokes_to_remove(i)).primitive_type = -2;
    end                      
end

function line = merge(vec1, vec2, w1, w2)
    norm = w1 + w2;
    w1 = sqrt(w1/norm);
    w2 = sqrt(w2/norm);
    
    x = [vec1([1,2])'; vec2([1,2])'];
    y = [vec1([3,4])'; vec2([3,4])'];


    [~, ind] = max([max(x)-min(x) max(y)-min(y)]);
    if (ind == 1)
        [f, ~] = fit(x,y,'poly1');
    else
        [f, ~] = fit(y,x,'poly1');
    end
      options = optimoptions('lsqnonlin', 'Algorithm', 'trust-region-reflective',...
                            'Display', 'off');  
                            lb =[];
    ub =[];
                        
    if (ind == 1)
        [x1, ind_x1] = min(x);
        [x2, ind_x2] = max(x);
        y1 = f(x1);
        y2 = f(x2);  

        
        
        param(1) = (y2-y1)/(x2 - x1);
        param(2) = -param(1)*x1 + y1;
        
%         plot([x1:x2], param(1)*[x1:x2] + param(2), 'r:')
        
        [param] = lsqnonlin(@findBestWeightedLineY, param, lb, ub, options);
        
%         plot([x1:x2], param(1)*[x1:x2] + param(2), 'r')
        
        y1 = param(1)*x1 + param(2);
        y2 = param(1)*x2 + param(2);
        
        % Ensure that the line endpoints are in the same order as in the
        % polyline:
        if (ind_x1 < ind_x2)
            line = [x1 x2 y1 y2];
        else
            line = [x2 x1 y2 y1];
        end
        
    else
        [y1, ind_y1] = min(y);
        [y2, ind_y2] = max(y);
        x1 = f(y1);
        x2 = f(y2);
        
        param(1) = (x2 - x1)/(y2-y1);
        param(2) = -param(1)*x1 + y1;
%          plot(param(1)*[y1:y2] + param(2), [y1:y2], 'r:')
        [param] = lsqnonlin(@findBestWeightedLineX, param, lb, ub, options);
%          plot(param(1)*[y1:y2] + param(2), [y1:y2], 'r')
        x1 = param(1)*y1 + param(2);
        x2 = param(1)*y2 + param(2);
        
        % Ensure that the line endpoints are in the same order as in the
        % polyline:
        if (ind_y1 < ind_y2)
            line = [x1 x2 y1 y2];
        else
            line = [x2 x1 y2 y1];
        end
    end



    function distance = findBestWeightedLineY(param)
        distance(1) = (vec1(3) - (param(1)*vec1(1) + param(2)))*w1;
        distance(2) = (vec1(4) - (param(1)*vec1(2) + param(2)))*w1;
        distance(3) = (vec2(3) - (param(1)*vec2(1) + param(2)))*w2;
        distance(4) = (vec2(4) - (param(1)*vec2(2) + param(2)))*w2;
    end

    function distance = findBestWeightedLineX(param)
        distance(1) = (vec1(1) - (param(1)*vec1(3) + param(2)))*w1;
        distance(2) = (vec1(2) - (param(1)*vec1(4) + param(2)))*w1;
        distance(3) = (vec2(1) - (param(1)*vec2(3) + param(2)))*w2;
        distance(4) = (vec2(2) - (param(1)*vec2(4) + param(2)))*w2;
    end
%     figure(12);        
%     hold off;
%     plot(vec1([1,2]), vec1([3,4]), 'b');        
%     hold on;
%     plot(vec2([1,2]), vec2([3,4]), 'g');
% 
%     axis equal
%     plot(line([1,2]), line([3,4]), 'r');
%       

end





function [begin, end_coordinate,vec1,vec2,dist] = findEndpointsCommonRegion(vec1, vec2)
    


   % 1. Align the strokes

    dir1 = vec1([2,4]) - vec1([1,3]);
    dir2 = vec2([2,4]) - vec2([1,3]);
    if dot(dir1./norm(dir1), dir2./norm(dir2)) < 0
        temp = vec2([2,4]);
        vec2([2,4]) = vec2([1,3]);
        vec2([1,3]) = temp;
    end

    % 2. Find the shortes distances at each of endpoints
    % Beginning:
    d1(1) = findLineSegmentPointDistance(vec2, vec1([1,3]));

    d1(2) = findLineSegmentPointDistance(vec1, vec2([1,3]));

    [d1_,ind] = min(d1);
    begins = [vec1([1,3]); vec2([1,3])];
    begin = begins(ind,:);
    
   % end:
    d1(1) = findLineSegmentPointDistance(vec2,vec1([2,4]));

    d1(2) = findLineSegmentPointDistance(vec1,vec2([2,4]));

    [d2_,ind] = min(d1);
    end_poitns = [vec1([2,4]); vec2([2,4])];
    end_coordinate = end_poitns(ind,:);
    
    dist = (d1_ + d2_)*0.5;

        
end


