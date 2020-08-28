clear all;
folder = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\curves17052020\results\results_medium';
folder_images = 'C:\Users\yulia\Research\Data\sketches_json_first_viewpoint';
folder_files = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\_400';

global folder_save;
global sketch_height;
sketch_height = 512;
global folder_save_base;
folder_save_base = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\curves17052020\results\medium_frames_rough';
global designer;
global object;
files = dir(folder);
files = files(~[files(:).isdir]);
files = {files.name};
% folder_to = 'D:\Documents\ForVideo\OriginalView';

for i = 1:length(files)
   close all
   
   filename = files{i}
   filename_ = strrep(filename, '.json', '');
   members = split(filename_,'_');
   designer = members{1};
   
   object_parts = members(2:(end-2));
   object = object_parts{1};
   for j = 2:length(object_parts)
      object = sprintf('%s_%s', object, object_parts{j});
   end
   
   view = 'view1';
%    
%    filepath_sketch_img = fullfile(folder_images,...
%                                    designer,...
%                                    object,...
%                                    [view '_' 'concept' '_opaque.png']);
%     
%    img = readSketchImg(filepath_sketch_img, true);                           
    
   try
       [ strokes_topology_curves, intersections, cam_param] = ...
            readReconstructionJsonFelix( fullfile(folder, filename));   
   catch
       continue;
   end
      
   [ strokes_topology_rough, ~, ~] = ...
            readReconstructionJson( fullfile(folder_files,designer, object, view, filename)); 
        
    if length(strokes_topology_rough) ~= length(strokes_topology_curves)
       continue; 
    end
        
    for ii = 1:length(strokes_topology_rough)
        
       
        if (strokes_topology_curves(ii).primitive_type ~= 0) & (strokes_topology_curves(ii).primitive_type ~= -2)
            continue;
        end
        disp(ii)
        strokes_topology_curves(ii).points3D = strokes_topology_rough(ii).points3D;
        if ~isnan(strokes_topology_curves(ii).points3D)
            strokes_topology_curves(ii).primitive_type = 0;
            strokes_topology_curves(ii).depth_assigned = true;
        end
    end
    
    
%     file_from = fullfile(folder_files, designer, object, view, 'sketch.svg');
%     file_to = fullfile(folder_to, [designer '_' object '.svg']);
%     copyfile(file_from,file_to );
%     copyfile()
%         
   [strokes_topology_curves, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology_curves, intersections, cam_param, []);
    
   folder_save = fullfile(folder_save_base, filename_);
   
   if ~exist(folder_save, 'dir')
       mkdir(folder_save);
   end
   
%    objectLFSquareAndSaveSVGFrames(strokes_topology_curves, cam_param, 'LFS')
    try
        rotate3DAndSaveSVGFrames(strokes_topology_curves, cam_param, '360');
    catch e
        disp(e);
    end
%    folder_svg = fullfile(folder_save, 'animation', 'LFS');

%     if ~exist(folder_svg, 'dir')
%         rmdir(folder_svg);
%     end

%    objectLFEffectAndSaveSVGFrames(strokes_topology_curves, cam_param, 'lf')
%    objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology_curves, cam_param, 'scaleXYZ')    
end

% [ strokes_topology_curves, intersections, cam_param] = ...
% readReconstructionJsonFelix( fullfile(folder, 'Professional2_vacuum_cleaner_bestScore_full.json' ));
% [strokes_topology_curves, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology_curves, intersections, cam_param, []);
% folder_save = fullfile(folder_save_base, 'Professional2_vacuum_cleaner');
% 
% objectLFEffectAndSaveSVGFrames(strokes_topology_curves, cam_param, 'lf')
% objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology_curves, cam_param, 'scaleXYZ')
% % rotate3DAndSaveSVGFrames(strokes_topology_curves, cam_param, '360');
% 
% [ strokes_topology_curves, intersections, cam_param] = ...
% readReconstructionJsonFelix( fullfile(folder, 'Professional6_mouse_bestScore_full.json' ));
% [strokes_topology_curves, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology_curves, intersections, cam_param, []);
% folder_save = fullfile(folder_save_base, 'Professional6_mouse');
% % objectLFEffectAndSaveSVGFrames(strokes_topology_curves, cam_param, 'lf')
% objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology_curves, cam_param, 'scaleXYZ')
% % rotate3DAndSaveSVGFrames(strokes_topology_curves, cam_param, '360');
% 
% [ strokes_topology_curves, intersections, cam_param] = ...
% readReconstructionJsonFelix( fullfile(folder, 'student8_house_bestScore_full.json' ));
% [strokes_topology_curves, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology_curves, intersections, cam_param, []);
% folder_save = fullfile(folder_save_base, 'student8_house');
% %objectLFEffectAndSaveSVGFrames(strokes_topology_curves, cam_param, 'lf')
% objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology_curves, cam_param, 'scaleXYZ')
% % rotate3DAndSaveSVGFrames(strokes_topology_curves, cam_param, '360');