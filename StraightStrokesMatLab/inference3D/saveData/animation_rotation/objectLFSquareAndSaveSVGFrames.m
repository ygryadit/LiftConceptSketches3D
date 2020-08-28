function objectLFSquareAndSaveSVGFrames(strokes_topology, cam_param, folder_name)

% step = 30;
num_views = 10;
        
fl = cam_param.K(1,1);
% values_fl = [(fl-num_views*step):step:(fl+num_views*step)];

% step_c = 0.04;
c = cam_param.C;
cam_param.view_dir = cam_param.view_dir./norm(cam_param.view_dir);

vec = cross(cam_param.view_dir, [0 0 1]');
up = cross(cam_param.view_dir,vec);
vec = cross(cam_param.view_dir, up);

step_c = 0.005;
c = cam_param.C;

tt = num_views*2 +1;
mm = num_views+1;



alpha = lsqnonlin(@distCenter, norm(cam_param.view_dir+cam_param.C));

focal_point = cam_param.C + alpha * cam_param.view_dir;

l = norm(c+focal_point);
step_c = l*0.001;


steps_h = num_views*step_c:-step_c:-num_views*step_c;
steps_w = num_views*step_c:-step_c:-num_views*step_c;

l = 0;
for i = 1:(num_views*2 + 1)
    for j = 1:(num_views*2 + 1)
        values(1, i, j) = cam_param.C(1) + vec(1)*steps_h(i) + up(1)*steps_w(j);
        values(2, i, j) = cam_param.C(2) + vec(2)*steps_h(i) + up(2)*steps_w(j);
        values(3, i, j) = cam_param.C(3) + vec(3)*steps_h(i) + up(3)*steps_w(j);
        c = values(:,i,j);
%         figure(4);
%         plot3(c(1),c(2),c(3),'*');
        moveCamAndSaveSVG(strokes_topology, cam_param, c, focal_point, l, folder_name);    
        l = l+1;
    end
end
    


%% 
function f = distCenter(alpha)
    f = cam_param.C + cam_param.view_dir*alpha';
end

end





function points = makeCircle(num_points, center, radius)
    s1 = num_points; 
    t1 = 0:(pi/(5*s1)):2*pi; % s1 is the number of points I want 
    
    r1 = radius; 
    x1 = center(1);
    y1 = center(2);
    x1unit = r1 * cos(t1) + x1; 
    y1unit = r1 * sin(t1) + y1;
    
    points(1,:)=x1unit(1:10:s1*10); %points around on circle
    points(2,:)=y1unit(1:10:s1*10);

end


function moveCamAndSaveSVG(strokes_topology, cam_param, c, focal_point, frame_num, folder_name)
global folder_save;
 save_color = false;

folder_svg = fullfile(folder_save, 'animation', folder_name);

if ~exist(folder_svg, 'dir')
    mkdir(folder_svg);
end

folder_black = fullfile(folder_svg, 'black');
if ~exist(folder_black , 'dir')
    mkdir(folder_black );
end

if save_color 
    folder_color = fullfile(folder_svg, 'color');
    if ~exist(folder_color, 'dir')
        mkdir(folder_color);
    end
end

file_svg = fullfile(folder_black, sprintf('sketch_black%03d.svg', frame_num))
if save_color 
    file_svg_color = fullfile(folder_color, sprintf('sketch_color%03d.svg', frame_num))
end

fidSVG = fopen(file_svg, 'w');
fidSVG = svgHead(fidSVG);

if save_color 
    fidSVGColor = fopen(file_svg_color, 'w');
    fidSVGColor = svgHead(fidSVGColor);
end

color_black = [0 0 0];

% cam_param.K(1,1) = fl;
% cam_param.K(2,2) = fl;
% cam_param.C = c;
% cam_param.P =  cam_param.K *[ cam_param.R -cam_param.R*cam_param.C];

% global filepath_sketch_img;
% img = readSketchImg(filepath_sketch_img, true);
img = ones(512,512,3);

up =[0 0 1];
cam_param.C = c;
cam_pos = cam_param.C;
cam_pos = reshape(cam_pos, 1, 3);

% cam_param.R
cam_param.R = rotationMatrixFromView(cam_pos, focal_point', up);
% cam_param.R

cam_param.P =  cam_param.K *[ cam_param.R -cam_param.R*cam_param.C];

%% Rotate and write
inds_assigned_depth = find(cat(1,strokes_topology(:).depth_assigned));
for i = 1:length(inds_assigned_depth)
    strk_num = inds_assigned_depth(i);
    
    points3D = strokes_topology(strk_num).points3D(1,:)';
%     points3D(1) = points3D(1)*sc;
%     strokes_topology(strk_num).points3D(1,:) = points3D';
    point2D = cam_param.P*[points3D; 1.0];
    point2D = point2D(1:2)./point2D(3);
    
    fidSVG = addPointStrokeBeginSVG(fidSVG, point2D);
    if save_color 
        fidSVGColor = addPointStrokeBeginSVG(fidSVGColor, point2D);
    end
    
    for ip = 2:size(strokes_topology(strk_num).points3D,1)
        points3D = strokes_topology(strk_num).points3D(ip,:)';
%         points3D(1) = points3D(1)*sc;
%         strokes_topology(strk_num).points3D(ip,:) = points3D';
        
        point2D = cam_param.P*[points3D; 1.0];
        point2D = point2D(1:2)./point2D(3);
        fidSVG = addPointSVG(fidSVG, point2D);
        if save_color 
            fidSVGColor = addPointSVG(fidSVGColor, point2D);
        end
    end   
    
    fidSVG = strokeEnd(fidSVG, color_black, strokes_topology(strk_num).mean_pressure);
    
    if save_color 
        fidSVGColor = strokeEnd(fidSVGColor, selectColor(strokes_topology(strk_num).line_group), strokes_topology(strk_num).mean_pressure);
    end
end

% reproject3Dto2D(img, cam_param, strokes_topology,[],2);


svgTail(fidSVG);
if save_color 
    svgTail(fidSVGColor);
end

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