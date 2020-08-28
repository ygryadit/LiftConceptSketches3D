function rotate3DAndSaveSVGFrames(strokes_topology, cam_param, folder_name)
    step = 2.5;
    num_views = 360/step;
    angle = 0;

    
    for i = 1:num_views
        rotate3DAndSaveSVG(strokes_topology, cam_param, angle, i-1, folder_name);
        angle = angle + step;
    end
end


function rotate3DAndSaveSVG(strokes_topology, cam_param, angle, frame_num, folder_name)
global folder_save;

folder_svg = fullfile(folder_save, 'animation', folder_name);

if ~exist(folder_svg, 'dir')
    mkdir(folder_svg);
end



folder_black = fullfile(folder_svg, 'black');
if ~exist(folder_black , 'dir')
    mkdir(folder_black );
end

folder_color = fullfile(folder_svg, 'color');
if ~exist(folder_color, 'dir')
    mkdir(folder_color);
end

file_svg = fullfile(folder_black, sprintf('sketch_black%03d.svg', frame_num));
if frame_num == 0
    global folder_save_base;
    global designer;
    global object;
    file_svg_ref_view = fullfile([folder_save_base '_ref_view'], [designer '_' object '.svg']);
    
    if ~exist([folder_save_base '_ref_view'])
        mkdir([folder_save_base '_ref_view']);
    end
    
    fidSVGREFView = fopen(file_svg_ref_view, 'w');
    fidSVGREFView = svgHead(fidSVGREFView);
end


fidSVG = fopen(file_svg, 'w');
fidSVG = svgHead(fidSVG);



color_black = [0 0 0];

R = rotz_mine(angle);
%% Rotate and write
inds_assigned_depth = find(cat(1,strokes_topology(:).depth_assigned));
for i = 1:length(inds_assigned_depth)
    strk_num = inds_assigned_depth(i);
    
    point2D = getProjection(strokes_topology(strk_num).points3D(1,:), R, cam_param.P);
    
    fidSVG = addPointStrokeBeginSVG(fidSVG, point2D);
    if frame_num == 0
        fidSVGREFView = addPointStrokeBeginSVG(fidSVGREFView, point2D);
    end
    
    for ip = 2:size(strokes_topology(strk_num).points3D,1)
        point2D = getProjection(strokes_topology(strk_num).points3D(ip,:), R, cam_param.P);
        fidSVG = addPointSVG(fidSVG, point2D);
        if frame_num == 0
            fidSVGREFView = addPointSVG(fidSVGREFView, point2D);
        end

    end   
    
    fidSVG = strokeEnd(fidSVG, color_black, strokes_topology(strk_num).mean_pressure);
    if frame_num == 0
        fidSVGREFView = strokeEnd(fidSVGREFView, color_black, strokes_topology(strk_num).mean_pressure);
    end
    
%     fidSVGColor = strokeEnd(fidSVGColor, selectColor(strokes_topology(strk_num).line_group), strokes_topology(strk_num).mean_pressure);
    
end

if frame_num == 0
    svgTail(fidSVGREFView);
end
svgTail(fidSVG);

end

function color = selectColor(line_group)
   if line_group == 1 
        color = uint8([1.0 0 0.5]*255);
   elseif line_group == 2
        color = uint8([0 0.8 0.4]*255);   
   elseif line_group == 3 
        color = uint8([0.4 0.3 1.0]*255);
   else
        color = uint8([0.8 0.8 0.0]*255);
   end
end

function point2D = getProjection(point3D, R, P)
    point3D = [R*point3D'; 1.0];

    point2D = P*point3D;
    point2D = point2D(1:2)./point2D(3);
end

function fidSVG = addPointStrokeBeginSVG(fidSVG, point2D)     
   str = sprintf('<path d="M %.5f %.5f \n',...
                   point2D(1), ...
                   point2D(2));   
   fwrite(fidSVG, str);
end

function fidSVG = addPointSVG(fidSVG, point2D)
  str = sprintf('L %.5f %.5f \n', ...
               point2D(1), ...
               point2D(2));   
  fwrite(fidSVG, str);
end

function fidSVG = strokeEnd(fidSVG, color, avg_pressure)
   pen_width = 1.5;
   str = sprintf('" fill="none" stroke="rgb(%d, %d, %d)" stroke-width="%.5f" stroke-opacity="%.3f" /> \n',...
                color(1),...
                color(2),...
                color(3),...
                pen_width*avg_pressure,...
                avg_pressure);   
            
    fwrite(fidSVG, str);
end

function fidSVG = svgHead(fidSVG)
global sketch_height;
str = sprintf('<?xml version="1.0" encoding="utf-8" ?>\n'); 
fwrite(fidSVG, str);
str = sprintf('<svg baseProfile="full" height="%d" version="1.1" viewBox="0,0,%d,%d" width="%d" xmlns="http://www.w3.org/2000/svg" xmlns:ev="http://www.w3.org/2001/xml-events" xmlns:xlink="http://www.w3.org/1999/xlink"><defs><style type="text/css"><![CDATA[\n', ...
    sketch_height,sketch_height,sketch_height,sketch_height);
fwrite(fidSVG, str);
str = sprintf('\t.background { fill: white; }\n');
fwrite(fidSVG, str);
str = sprintf('\t.line { stroke: firebrick; stroke-width: .1mm; }\n');
fwrite(fidSVG, str);
str = sprintf('\t.blacksquare { fill: indigo; }\n');
fwrite(fidSVG, str);
str = sprintf('\t.whitesquare { fill: white; }\n');
fwrite(fidSVG, str);
str = sprintf(']]></style></defs>\n');
fwrite(fidSVG, str);
end


function svgTail(fidSVG)
str = sprintf('</svg>');
fwrite(fidSVG, str);

fclose(fidSVG);

end