function inds_corners = getIndsCorners(inds_vp,...
                            inds_corners,...
                            vpgr1,...
                            vpgr2,...
                            inds_strks_1,...
                            inds_strks_2,...
                            strokes_topology,...
                            intersections)
% global filepath_sketch_img;
% img = readSketchImg(filepath_sketch_img);
% figure(2);
% imshow(img);
% hold on;

for i = inds_vp
    if ismember(i, inds_corners)
        continue;
    end
       
    
    
    
    
    indcs_intrsctns_1 = strokes_topology(inds_strks_1(i)).indcs_intrsctns;
    mask = ~intersections.collinear((indcs_intrsctns_1));    
    inds_intrsctng_strks_1 = setdiff(strokes_topology(inds_strks_1(i)).indcs_intrsctng_strks(mask), inds_strks_2(i));    
    
    indcs_intrsctns_2 = strokes_topology(inds_strks_2(i)).indcs_intrsctns;
    mask = ~intersections.collinear((indcs_intrsctns_2));    
    inds_intrsctng_strks_2 = setdiff(strokes_topology(inds_strks_2(i)).indcs_intrsctng_strks(mask), inds_strks_1(i));        
    clear('mask');
    
    
    inds_intrsctng_strks_1 = exculudeCollinear(inds_strks_2(i), inds_intrsctng_strks_1, inds_strks_1, inds_strks_2);
    inds_intrsctng_strks_2 = exculudeCollinear(inds_strks_1(i), inds_intrsctng_strks_2, inds_strks_1, inds_strks_2);
    
   
    
%     plot([strokes_topology(inds_strks_1(i)).points2D(:).x], ...
%          [strokes_topology(inds_strks_1(i)).points2D(:).y], 'm');
%     plot([strokes_topology(inds_strks_2(i)).points2D(:).x], ...
%          [strokes_topology(inds_strks_2(i)).points2D(:).y], 'k');
     
    
    mask_1 = [strokes_topology(inds_intrsctng_strks_1).line_group] == vpgr2;
    mask_2 = [strokes_topology(inds_intrsctng_strks_2).line_group] == vpgr1;
    
    inds_intrsctng_strks_1 = inds_intrsctng_strks_1(mask_1);
    inds_intrsctng_strks_2 = inds_intrsctng_strks_2(mask_2);
    [p,q] = meshgrid(inds_intrsctng_strks_1, inds_intrsctng_strks_2);
    pairs = [p(:) q(:)];

%     for ii = 1:length(inds_intrsctng_strks_1)
%         plot([strokes_topology(inds_intrsctng_strks_1(ii)).points2D(:).x], ...
%             [strokes_topology(inds_intrsctng_strks_1(ii)).points2D(:).y], 'b:');
%     end
%     
%     for ii = 1:length(inds_intrsctng_strks_2)
%          plot([strokes_topology(inds_intrsctng_strks_2(ii)).points2D(:).x], ...
%              [strokes_topology(inds_intrsctng_strks_2(ii)).points2D(:).y], 'r:');
%     end
    
    % Points of intersectiosn of two additional strokes:
    mask_cycle_1 = ismember(inds_strks_1, pairs(:,1)) & ismember(inds_strks_2, pairs(:,2));
    mask_cycle_2 = ismember(inds_strks_2, pairs(:,1)) & ismember(inds_strks_1, pairs(:,2));
    mask_cycle = mask_cycle_1 | mask_cycle_2;
    inds_cycle = find(mask_cycle);
    
    
    strks_pairs = intersections.strokes_indices(inds_cycle,:);
    %rearrange so that the first stoke is towards second vp, and the second
    %towards first:
    inds = [strokes_topology(strks_pairs(:,1)).line_group] ~= vpgr2;
    temp = strks_pairs(inds,1);
    strks_pairs(inds,1) = strks_pairs(inds,2);
    strks_pairs(inds,2) = temp;
    
    
    [~,locb] = ismember(strks_pairs(:,1), ...
                        strokes_topology(inds_strks_1(i)).indcs_intrsctng_strks);
    
    i1 = strokes_topology(inds_strks_1(i)).indcs_intrsctns( locb);
    
    [~,locb] = ismember(strks_pairs(:,2), ...
                        strokes_topology(inds_strks_2(i)).indcs_intrsctng_strks);
                    
    
    i2 = strokes_topology(inds_strks_2(i)).indcs_intrsctns(locb);
                    
    
    
    
    if ~isempty(inds_cycle)
        inds_corners_ = [ repmat(i, length(inds_cycle), 1)...
                            i1...
                            inds_cycle ...
                            i2];
                    
%         plot(intersections.coordinates2D(inds_corners_(:), 1), ...
%              intersections.coordinates2D(inds_corners_(:), 2), '*');
        inds_corners = [inds_corners; inds_corners_];
    end
    
end


function list2 = exculudeCollinear(i, list2, inds_strks_1, inds_strks_2)
    %remove collinear strokes
    list_1 = repmat(i,length(list2),1);
    
    mask1 = ismember(inds_strks_1,list_1) & ...
           ismember(inds_strks_2,list2);
    
    strks1 = inds_strks_2(mask1);
       
       
    mask2 = ismember(inds_strks_2,list_1) & ...
            ismember(inds_strks_1,list2);
                   ismember(inds_strks_2,list2);
    
    strks2 = inds_strks_1(mask2);
    
    
    
    list2 = setdiff(list2, [strks1; strks2]);
 
end

end


% function mask_collinear = exculudeCollinear(i, list2, inds_strks_1, inds_strks_2)
%     %remove collinear strokes
%     list_1 = repmat(i,length(list2),1);
%     
%     mask1 = ismember(inds_strks_1,list_1) & ...
%            ismember(inds_strks_2,list2);
%     ind1 = find(mask1)    ;
%     if ~isempty(ind1)
%     for ii = ind1'
%          plot([strokes_topology(inds_strks_2(ii)).points2D(:).x], ...
%              [strokes_topology(inds_strks_2(ii)).points2D(:).y], 'k');
%     end
%     end
%        
%     mask2 = ismember(inds_strks_2,list_1) & ...
%             ismember(inds_strks_1,list2);
%                    ismember(inds_strks_2,list2);
%     
%     ind1 = find(mask2);
%     if ~isempty(ind1)
%     for ii = ind1'
%          plot([strokes_topology(inds_strks_2(ii)).points2D(:).x], ...
%              [strokes_topology(inds_strks_2(ii)).points2D(:).y], 'k');
%     end
%     end
%     
%     mask_collinear = mask1 | mask2;
%     
%  
% end