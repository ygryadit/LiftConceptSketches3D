function [vps_selected,...
    inds_axis,...
    inds_inactive_lines] = computeRotatedSetsOfParallelLines(vp,p,All_lines,w,h,img)

    
    vp =reshape(vp', 2, 3)';
    %% Find remainimg lines:

    [~, linemem] = max(p,[],2);            
    grp4=linemem==4;
    
    lines_active = All_lines(grp4, :);
    inds_lns = find(grp4);
    
    
    vps_selected =[];
    inds_axis = [];
    inds_inactive_lines ={};
    i = 1;
    while (~isempty(lines_active))
        try
        [vps_selected_,...
         inds_axis_,...
         inds_inactive_lines_,...
         inds_active_lines] = computeOneRotatedSetOfParallelLines(vp,lines_active,w,h,img);
        catch e
           disp(''); 
           rethrow(e);
        end
        if length(inds_inactive_lines_) > 4
            vps_selected(i,:) = vps_selected_;
            inds_axis(i) = inds_axis_;
            inds_inactive_lines{i} = inds_lns(inds_inactive_lines_);
        else
            break;
        end
     
        i = i+1;
        lines_active = lines_active(inds_active_lines,:);
        inds_lns = inds_lns(inds_active_lines);
    end
end

