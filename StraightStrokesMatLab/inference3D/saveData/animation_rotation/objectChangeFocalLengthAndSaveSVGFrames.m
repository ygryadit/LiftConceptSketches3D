function objectChangeFocalLengthAndSaveSVGFrames(strokes_topology, cam_param, folder_name)

step = 30;
num_views = 10;
        
fl = cam_param.K(1,1);
values_fl = [(fl-num_views*step):step:(fl+num_views*step)];

step_c = 0.04;
c = cam_param.C;
vdir = cam_param.view_dir;
if -vdir(1)*step_c ~= 0
    values(1,:) = [(c(1)+vdir(1)*num_views*step_c):-vdir(1)*step_c:(c(1)-vdir(1)*num_views*step_c)];
else
    values(1,:) = ones(1,21)*c(1);   
end

if -vdir(2)*step_c ~= 0
    values(2,:) = [(c(2)+vdir(2)*num_views*step_c):-vdir(2)*step_c:(c(2)-vdir(2)*num_views*step_c)];
else
    values(2,:) = ones(1,21)*c(2);  
end
if -vdir(3)*step_c ~= 0
    values(3,:) = [(c(3)+vdir(3)*num_views*step_c):-vdir(3)*step_c:(c(3)-vdir(3)*num_views*step_c)];
else
    values(3,:) = ones(1,21)*c(3);
end


values_fl = values_fl([11:21 20:-1:1 2:11]);
values = values(:, [11:21 20:-1:1 2:11]);

% values_fl = values_fl([11:21 20:-1:11]);
% values = values(:, [11:21 20:-1:11]);

for i = 1:length(values_fl)
    fl = values_fl(i);
    c = values(:,i);
    changeFocalLengthAndSaveSVG(strokes_topology, cam_param, fl, c, i-1, folder_name);    
end

end


function changeFocalLengthAndSaveSVG(strokes_topology, cam_param, fl,c, frame_num, folder_name)
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

file_svg = fullfile(folder_black, sprintf('sketch_black%03d.svg', frame_num))
file_svg_color = fullfile(folder_color, sprintf('sketch_color%03d.svg', frame_num))

fidSVG = fopen(file_svg, 'w');
fidSVG = svgHead(fidSVG);

fidSVGColor = fopen(file_svg_color, 'w');
fidSVGColor = svgHead(fidSVGColor);

color_black = [0 0 0];

cam_param.K(1,1) = fl;
cam_param.K(2,2) = fl;
cam_param.C = c;
cam_param.P =  cam_param.K *[ cam_param.R -cam_param.R*cam_param.C];

% global filepath_sketch_img;
% img = readSketchImg(filepath_sketch_img, true);
img = ones(512,512,3);
reproject3Dto2D(img, cam_param, strokes_topology,[],2);

%% Rotate and write
inds_assigned_depth = find(cat(1,strokes_topology(:).depth_assigned));
for i = 1:length(inds_assigned_depth)
    strk_num = inds_assigned_depth(i);
    
    point2D = cam_param.P*[strokes_topology(strk_num).points3D(1,:)'; 1.0];
    point2D = point2D(1:2)./point2D(3);
    
    fidSVG = addPointStrokeBeginSVG(fidSVG, point2D);
    fidSVGColor = addPointStrokeBeginSVG(fidSVGColor, point2D);
    
    for ip = 2:size(strokes_topology(strk_num).points3D,1)
        point2D = cam_param.P*[strokes_topology(strk_num).points3D(ip,:)'; 1.0];
        point2D = point2D(1:2)./point2D(3);
        fidSVG = addPointSVG(fidSVG, point2D);
        fidSVGColor = addPointSVG(fidSVGColor, point2D);
    end   
    
    fidSVG = strokeEnd(fidSVG, color_black, strokes_topology(strk_num).mean_pressure);
    
    
    fidSVGColor = strokeEnd(fidSVGColor, selectColor(strokes_topology(strk_num).line_group), strokes_topology(strk_num).mean_pressure);
    
end

svgTail(fidSVG);
svgTail(fidSVGColor);

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