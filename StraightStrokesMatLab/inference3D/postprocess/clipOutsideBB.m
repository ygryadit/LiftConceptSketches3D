function strokes_topology = clipOutsideBB(strokes_topology, intersections, cam_param)
    
    [bound_box, intersections_active_coord3D, inds] = bbOFIntersections(intersections);

    inds_assigned_depth = find(cat(1,strokes_topology(:).depth_assigned) & ...
                               (cat(1,strokes_topology(:).primitive_type) == 0));
    img = ones(512,512,3);
   reproject3Dto2D(img, cam_param, strokes_topology,[],2);
   
   close(figure(4))
   strokes_assigned = strokes_topology(cat(1, strokes_topology(:).depth_assigned));
   plotStrokesTopology(4, strokes_assigned);
%    plot3(intersections_active_coord3D(:,1), ...
%        intersections_active_coord3D(:,2),...
%        intersections_active_coord3D(:,3), '*');
   
   plot3([bound_box(1,1),bound_box(1,2)],...
          [bound_box(2,1),bound_box(2,1)],...
          [bound_box(3,1),bound_box(3,1)], 'LineWidth', 2);
   plot3([bound_box(1,1),bound_box(1,2)],...
          [bound_box(2,2),bound_box(2,2)],...
          [bound_box(3,1),bound_box(3,1)], 'LineWidth', 2);   
   plot3([bound_box(1,1),bound_box(1,2)],...
          [bound_box(2,2),bound_box(2,2)],...
          [bound_box(3,2),bound_box(3,2)], 'LineWidth', 2);  
   plot3([bound_box(1,1),bound_box(1,2)],...
         [bound_box(2,1),bound_box(2,1)],...
         [bound_box(3,2),bound_box(3,2)], 'LineWidth', 2);      
 
   plot3([bound_box(1,1),bound_box(1,1)],...
          [bound_box(2,1),bound_box(2,2)],...
          [bound_box(3,1),bound_box(3,1)], 'LineWidth', 2);
   plot3([bound_box(1,2),bound_box(1,2)],...
          [bound_box(2,1),bound_box(2,2)],...
          [bound_box(3,1),bound_box(3,1)], 'LineWidth', 2);
   plot3([bound_box(1,1),bound_box(1,1)],...
          [bound_box(2,1),bound_box(2,2)],...
          [bound_box(3,2),bound_box(3,2)], 'LineWidth', 2);   
   plot3([bound_box(1,2),bound_box(1,2)],...
         [bound_box(2,1),bound_box(2,2)],...
         [bound_box(3,2),bound_box(3,2)], 'LineWidth', 2);      
     
     
%      text(intersections_active_coord3D(:,1), ...
%        intersections_active_coord3D(:,2),...
%        intersections_active_coord3D(:,3), ...
%         cellstr(num2str(inds')));    
%    
   debug = false;
   
    for j = 1:length(inds_assigned_depth)
        strk_num = inds_assigned_depth(j);
        for i = 2:size(strokes_topology(strk_num).points3D,1)
            
            point3D_rev = strokes_topology(strk_num).points3D(i-1,:);
            point3D = strokes_topology(strk_num).points3D(i,:);
            
               
     
            if (~checkPointInsideBB(bound_box, point3D) & checkPointInsideBB(bound_box, point3D_rev)) ...
                 | ...
               (checkPointAfterBB(bound_box, point3D)  & checkPointBeforeBB(bound_box, point3D_rev))
                 
                if debug
                    figure(4);
                    hold off;
                    strokes_assigned = strokes_topology(cat(1, strokes_topology(:).depth_assigned));
                    plotStrokesTopology(4, strokes_assigned);

                    plot3(point3D(1), point3D(2), point3D(3), '*');
                    plot3(point3D_rev(1), point3D_rev(2), point3D_rev(3), 'o');
                end
            
                strokes_topology(strk_num).points3D = ...
                    clipStroke(strokes_topology(strk_num).points3D, bound_box, i);
                
                break;
            end
            
            if (checkPointInsideBB(bound_box, point3D) & ~checkPointInsideBB(bound_box, point3D_rev)) ...
                  | ...
               (checkPointBeforeBB(bound_box, point3D)  & checkPointAfterBB(bound_box, point3D_rev))
                if debug
                    figure(4);
                    hold off;
                    strokes_assigned = strokes_topology(cat(1, strokes_topology(:).depth_assigned));
                    plotStrokesTopology(4, strokes_assigned);

                    plot3(point3D(1), point3D(2), point3D(3), '*');
                    plot3(point3D_rev(1), point3D_rev(2), point3D_rev(3), 'o');
                end
                
                num_points = size(strokes_topology(strk_num).points3D,1);
                
                 strokes_topology(strk_num).points3D = ...
                    clipStroke(strokes_topology(strk_num).points3D(end:-1:1,:), bound_box, num_points-(i-2));                
                 break;
            end
        end
    end
    for j = 1:length(inds_assigned_depth)
        strk_num = inds_assigned_depth(j);
        for i = 2:size(strokes_topology(strk_num).points3D,1)
            
            point3D_rev = strokes_topology(strk_num).points3D(i-1,:);
            point3D = strokes_topology(strk_num).points3D(i,:);
            
               
     
            if (~checkPointInsideBB(bound_box, point3D) & checkPointInsideBB(bound_box, point3D_rev)) ...
                 | ...
               (checkPointAfterBB(bound_box, point3D)  & checkPointBeforeBB(bound_box, point3D_rev))
                 
                 if debug
                    figure(4);
                    hold off;
                    strokes_assigned = strokes_topology(cat(1, strokes_topology(:).depth_assigned));
                    plotStrokesTopology(4, strokes_assigned);

                    plot3(point3D(1), point3D(2), point3D(3), '*');
                    plot3(point3D_rev(1), point3D_rev(2), point3D_rev(3), 'o');
                end
%             
                strokes_topology(strk_num).points3D = ...
                    clipStroke(strokes_topology(strk_num).points3D, bound_box, i);
%                 
                break;
            end
            
            if (checkPointInsideBB(bound_box, point3D) & ~checkPointInsideBB(bound_box, point3D_rev)) ...
                 | ...
               (checkPointBeforeBB(bound_box, point3D)  & checkPointAfterBB(bound_box, point3D_rev))
                if debug
                    figure(4);
                    hold off;
                    strokes_assigned = strokes_topology(cat(1, strokes_topology(:).depth_assigned));
                    plotStrokesTopology(4, strokes_assigned);

                    plot3(point3D(1), point3D(2), point3D(3), '*');
                    plot3(point3D_rev(1), point3D_rev(2), point3D_rev(3), 'o');
                end
                %                 
                num_points = size(strokes_topology(strk_num).points3D,1);
                
                 strokes_topology(strk_num).points3D = ...
                    clipStroke(strokes_topology(strk_num).points3D(end:-1:1,:), bound_box, num_points-(i-2));                
                 break;
            end
        end
    end
    
    

    reproject3Dto2D(img, cam_param, strokes_topology,[],3);
end


function isInside = checkPointBeforeBB(bound_box, point3D)    
    d = 0.005;
    isInside =  ((point3D(1) < (bound_box(1,1) - d) ) | ...
                 (point3D(2) < (bound_box(2,1)- d) ) | ...
                 (point3D(3) < (bound_box(3,1) - d) ) );
             
end

function isInside = checkPointAfterBB(bound_box, point3D)    
     d = 0.005;
    isInside =  ((point3D(1) > (bound_box(1,2)+d)  ) | ...
                 (point3D(2) > (bound_box(2,2)+d)  ) | ...
                 (point3D(3) > (bound_box(3,2)+d)  ) );
             
end


function isInside = checkPointInsideBB(bound_box, point3D)    
    
    isInside =  ((point3D(1) >= bound_box(1,1)) & ...
                 (point3D(1) <= bound_box(1,2)) & ...
                 (point3D(2) >= bound_box(2,1)) & ...
                 (point3D(2) <= bound_box(2,2)) & ...
                 (point3D(3) >= bound_box(3,1)) & ...
                 (point3D(3) <= bound_box(3,2)));
             
end


function points3D = clipStroke(points3D, bound_box, ind_segment)
    seg = points3D((ind_segment-1):ind_segment, :);

    dir = seg(2,:) - seg(1,:);
    dir = dir/norm(dir);
    
   point3D =  points3D(ind_segment, :);
   point3D_start =  points3D(ind_segment-1, :);
   f = 0.001;
    
    x1 = point3D(1) < bound_box(1,1);
    x2 = point3D(1) > bound_box(1,2);
    
   
    if (x1) | (x2)
        if x1
            x = bound_box(1,1);
        else
            x = bound_box(1,2);
        end
        
        [y,z] = intersectionPointGivenX(dir, x, point3D_start);
        
%         points3D = [points3D(1:ind_segment-1,:)];
        
%         plot3(points3D(:,1), points3D(:,2), points3D(:,3), 'r');
%         
%         plot3(point3D_start(1),point3D_start(2),point3D_start(3),'*');
        
%         plot3(point3D_start(1)+dir(1),point3D_start(2)+dir(2),point3D_start(3)+dir(3),'*');
        
        points3D = [points3D(1:(ind_segment-1),:); x,y,z;  [x,y,z] + dir*f];
%         figure(4);
%         plot3(x,y,z,'*');
        
        return;
    end
 
    
    y1 =(point3D(2) < bound_box(2,1));
    y2 = (point3D(2) > bound_box(2,2));
    
    if y1|y2
        if y1
            y = bound_box(2,1);
        else
            y = bound_box(2,2);            
        end
        
        [x,z] = intersectionPointGivenY(dir, y, point3D_start);
%         points3D = [points3D(1:ind_segment-1,:)];
%         plot3(point3D_start(:,1), point3D_start(:,2), point3D_start(:,3), 'r');
        
        points3D = [points3D(1:(ind_segment-1),:); x,y,z; [x,y,z] + dir*f];
%         figure(4);
%         plot3(x,y,z,'*');
        
        return;
    end
    
    z1 =  (point3D(3) < bound_box(3,1));
    z2 =  (point3D(3) > bound_box(3,2));
    
    if z1 | z2
        if z1
            z = bound_box(3,1);
        else
            z = bound_box(3,2);            
        end
        [x,y] = intersectionPointGivenZ(dir, z, point3D_start);    
%         points3D = [points3D(1:ind_segment-1,:)];
%         plot3(point3D_start(:,1), point3D_start(:,2), point3D_start(:,3), 'r');
        points3D = [points3D(1:(ind_segment-1),:); x,y,z;  [x,y,z] + dir*f];
%         figure(4);
%         plot3(x,y,z,'*');
        
        return;
    end
    
    
%     if (point3D(1) < bound_box(1,1))
%         x = bound_box(1,1);
%         [y,z] = intersectionPointGivenX(dir, x, pointstart);
%         
%         points3D = [points3D(1:ind_segment-1,:); [x,y,z] + dir*0.2];
%         return;
%     end
%  
%     if (point3D(2) > bound_box(2,1))
%         y = bound_box(2,1);
%         [x,z] = intersectionPointGivenY(dir, y, pointstart);
%         
%         points3D = [points3D(1:ind_segment-1,:); [x,y,z] + dir*0.2];
%         return;
%     end
%     
%     if (point3D(3) > bound_box(3,1))
%         z = bound_box(3,1);
%         [x,y] = intersectionPointGivenZ(dir, z, pointstart);        
%         points3D = [points3D(1:ind_segment-1,:); [x,y,z] + dir*0.2];
%         return;
%     end
    
    
    
end

function [y,z] = intersectionPointGivenX(dir, x, pointstart)
dir = dir/norm(dir);
 y = pointstart(2) + dir(2)/dir(1)*(x - pointstart(1));
 z = pointstart(3) + dir(3)/dir(1)*(x - pointstart(1)); 
end

function [x,z] = intersectionPointGivenY(dir, y, pointstart)
 dir = dir/norm(dir);
 x = pointstart(1) + dir(1)/dir(2)*(y - pointstart(2));
 z = pointstart(3) + dir(3)/dir(2)*(y - pointstart(2)); 
end

function [x,y] = intersectionPointGivenZ(dir, z, pointstart)
 dir = dir/norm(dir);
 x = pointstart(1) + dir(1)/dir(3)*(z - pointstart(3));
 y = pointstart(2) + dir(2)/dir(3)*(z - pointstart(3)); 
end