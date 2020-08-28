function plot2DCurStrokeIntersectingStokes(intersections, ...
                                            cur_stroke,...
                                            strokes_topology)

                                        
    h2 = figure(2); 
    hold off;
    set(h2, 'Name', 'Intersection 2D');
    global filepath_sketch_img;
    img = readSketchImg(filepath_sketch_img);
    
    imshow(img);
    hold on;
    for i = 1:cur_stroke.ind
        if strokes_topology(i).depth_assigned
            plot(cat(1,strokes_topology(i).points2D.x),cat(1,strokes_topology(i).points2D.y), 'LineWidth',1, 'Color', 'b', 'LineStyle', ':');
        else
            plot(cat(1,strokes_topology(i).points2D.x),cat(1,strokes_topology(i).points2D.y), 'LineWidth',1, 'Color', 'r', 'LineStyle', ':');
        end
    end
    plot(cat(1,strokes_topology(cur_stroke.ind).points2D.x),cat(1,strokes_topology(cur_stroke.ind).points2D.y), 'LineWidth',1, 'Color', 'r');
   
   
    ic2D = cat(1,intersections(cur_stroke.indcs_intrsctns).coordinates2D);
    mask_likely = cat(1,intersections(cur_stroke.indcs_intrsctns).likely);
    ic2D = ic2D(mask_likely,:);
    
    plot( ic2D(: ,1),...
          ic2D(:,2),...
               'g*');
%     plot( ic2D(~mask_likely ,1),...
%           ic2D(~mask_likely,2),...
%           'b*');
           
    if ~isempty(cur_stroke.inds_intrsctns_eval)
        ic2D = cat(1,intersections(cur_stroke.inds_intrsctns_eval).coordinates2D);
        plot( ic2D(:,1),...
              ic2D(:,2),...
               'ro', 'LineWidth', 2);
    end
    
    
    
    strokes_mult_hypothesis = cur_stroke.inds_intrsctng_strks_eval_mltpl_cnddts;
    prev_strokes_assigned = cur_stroke.inds_intrsctng_strks_eval_actv;
    
    for stroke=strokes_mult_hypothesis'       
       plot(cat(1,strokes_topology(stroke).points2D.x),cat(1,strokes_topology(stroke).points2D.y), 'g:', 'LineWidth',2);
    end
    
    for stroke=prev_strokes_assigned'       
       plot(cat(1,strokes_topology(stroke).points2D.x),cat(1,strokes_topology(stroke).points2D.y), 'b-', 'LineWidth',2);
    end
    
    
    
%     legend('intersections conf.','mult hypoth', 'assigned', 'cur stroke');
    a = cur_stroke.inds_intrsctns_eval(:); 
    b = num2str(a); 
    c = cellstr(b);
    text(ic2D(:,1),...
         ic2D(:,2),...
         c, 'FontSize',14);   
     
    
    non_active = setdiff(cur_stroke.indcs_intrsctns(mask_likely), cur_stroke.inds_intrsctns_eval(:)); 
    if ~isempty(non_active)
    b = num2str(non_active); 
    c = cellstr(b);
    ic2D_ = cat(1,intersections(non_active).coordinates2D);
    text(ic2D_(:,1),...
         ic2D_(:,2),...
         c, 'FontSize',14, 'Color', 'green');   
    end
    
    %     legend('intersections conf.','mult hypoth', 'assigned', 'cur stroke');
    a = cur_stroke.inds_intrsctns_eval(:); 
    b = num2str(a); 
    c = cellstr(b);
    text(ic2D(:,1),...
         ic2D(:,2),...
         c, 'FontSize',14);   
     
end