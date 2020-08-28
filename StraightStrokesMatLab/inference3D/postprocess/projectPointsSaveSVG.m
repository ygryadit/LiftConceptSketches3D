function strokes_topology = projectPointsSaveSVG(strokes_topology, scaleX, scaleY, scaleZ, R_o, P, file_svg, file_svg_color)
    global save_color;
    color_black = [0 0 0];
    global sketch_height;
    sketch_height = 1024;
    fidSVG = fopen(file_svg, 'w');
    fidSVG = svgHead(fidSVG);
    shift = [256; 256];
    if save_color 
        fidSVGColor = fopen(file_svg_color, 'w');
        fidSVGColor = svgHead(fidSVGColor);
    end
    
    inds_assigned_depth = find(cat(1,strokes_topology(:).depth_assigned));
    for i = 1:length(inds_assigned_depth)
        strk_num = inds_assigned_depth(i);

        point3D = scaleStrokeDim(strokes_topology(strk_num), 1, scaleX, scaleY, scaleZ, R_o);        
        strokes_topology(strk_num).points3D(1,:) = point3D;
        point2D = projectPoint(point3D, P);
        point2D = point2D + shift;
          
        fidSVG = addPointStrokeBeginSVG(fidSVG, point2D);
        if save_color 
            fidSVGColor = addPointStrokeBeginSVG(fidSVGColor, point2D);
        end

        for ip = 2:size(strokes_topology(strk_num).points3D,1)
            point3D = scaleStrokeDim(strokes_topology(strk_num), ip, scaleX, scaleY, scaleZ, R_o);  
            strokes_topology(strk_num).points3D(ip,:) = point3D;
            point2D = projectPoint(point3D, P);
            point2D = point2D + shift;
            
            fidSVG = addPointSVG(fidSVG, point2D);
            if save_color 
                fidSVGColor = addPointSVG(fidSVGColor, point2D);
            end
        end   

        fidSVG = strokeEnd(fidSVG, color_black, strokes_topology(strk_num).mean_pressure, i);

        if save_color 
            fidSVGColor = strokeEnd(fidSVGColor, selectColor(strokes_topology(strk_num).line_group), strokes_topology(strk_num).mean_pressure, i);
        end
    end
    
    svgTail(fidSVG);
    if save_color 
        svgTail(fidSVGColor);
    end
end


function point3D = scaleStrokeDim(stroke, ind_point, scx, scy, scz, R_o)
    point3D = stroke.points3D(ind_point,:)';
    point3D(1) = point3D(1)*scx;  
    point3D(2) = point3D(2)*scy;  
    point3D(3) = point3D(3)*scz;  
    point3D = R_o*point3D;    
end


function point2D = projectPoint(point3D, P)
    point2D = P*[point3D; 1.0];
    point2D = point2D(1:2)./point2D(3);
end

function fidSVG = addPointSVG(fidSVG, point2D)
  str = sprintf('L %.5f %.5f \n', ...
               point2D(1), ...
               point2D(2));   
  fwrite(fidSVG, str);
end

function fidSVG = strokeEnd(fidSVG, color, avg_pressure,id)
   pen_width = 1.5;
   str = sprintf('" fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" id ="%d"/> \n',...
                color(1),...
                color(2),...
                color(3),...
                pen_width*avg_pressure,...
                avg_pressure,...
                id);   
            
    fwrite(fidSVG, str);
end