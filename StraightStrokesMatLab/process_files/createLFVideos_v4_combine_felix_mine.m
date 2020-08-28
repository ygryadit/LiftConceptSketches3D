clear all;
folder = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\curves17052020\results\results_good';
folder_images = 'C:\Users\yulia\Research\Data\sketches_json_first_viewpoint';
folder_files = 'C:\Users\yulia\Research\DesignSketch3D\results\v50';

global folder_save;
global sketch_height;
sketch_height = 512;
foler_save_base = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\curves17052020\results\good_frames_rough';

files = dir(folder);
files = files(~[files(:).isdir]);
files = {files.name};
% folder_to = 'D:\Documents\ForVideo\OriginalView';

for i = 1:length(files)
    close all
   filename = files{i};
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
       [ strokes_topology, intersections, cam_param] = ...
            readReconstructionJsonFelix( fullfile(folder, filename));   
   catch
       continue;
   end
   
     
   
%    launchOpenSketch_cluster_a2_400('C:\Users\yulia\Research\Data\sketches_json_first_viewpoint', designer, object, 'C:\Users\yulia\Research\DesignSketch3D\results\v50')
   [ strokes_topology_, ~, ~] = ...
            readReconstructionJson( fullfile(folder_files,designer, object, view, filename)); 
        
    for ii = 1:length(strokes_topology)
        strokes_topology(ii).primitive_type 
       
        if (strokes_topology(ii).primitive_type ~= 0) & (strokes_topology(ii).primitive_type ~= -2)
            continue;
        end
        disp(ii)
        strokes_topology(ii).points3D = strokes_topology_(ii).points3D;
        if ~isnan(strokes_topology(ii).points3D)
            strokes_topology(ii).primitive_type = 0;
            strokes_topology(ii).depth_assigned = true;
        end
    end
    
    
%     file_from = fullfile(folder_files, designer, object, view, 'sketch.svg');
%     file_to = fullfile(folder_to, [designer '_' object '.svg']);
%     copyfile(file_from,file_to );
%     copyfile()
%         
   [strokes_topology, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology, intersections, cam_param, []);
    
   folder_save = fullfile(foler_save_base, filename_);
   
   if ~exist(folder_save, 'dir')
       mkdir(folder_save);
   end
   
%    objectLFSquareAndSaveSVGFrames(strokes_topology, cam_param, 'LFS')
%    rotate3DAndSaveSVGFrames(strokes_topology, cam_param, '360');
   
%    folder_svg = fullfile(folder_save, 'animation', 'LFS');

%     if ~exist(folder_svg, 'dir')
%         rmdir(folder_svg);
%     end

   objectLFEffectAndSaveSVGFrames(strokes_topology, cam_param, 'lf')
%    objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology, cam_param, 'scaleXYZ')    
end

% [ strokes_topology, intersections, cam_param] = ...
% readReconstructionJsonFelix( fullfile(folder, 'Professional2_vacuum_cleaner_bestScore_full.json' ));
% [strokes_topology, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology, intersections, cam_param, []);
% folder_save = fullfile(foler_save_base, 'Professional2_vacuum_cleaner');
% 
% objectLFEffectAndSaveSVGFrames(strokes_topology, cam_param, 'lf')
% objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology, cam_param, 'scaleXYZ')
% % rotate3DAndSaveSVGFrames(strokes_topology, cam_param, '360');
% 
% [ strokes_topology, intersections, cam_param] = ...
% readReconstructionJsonFelix( fullfile(folder, 'Professional6_mouse_bestScore_full.json' ));
% [strokes_topology, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology, intersections, cam_param, []);
% folder_save = fullfile(foler_save_base, 'Professional6_mouse');
% % objectLFEffectAndSaveSVGFrames(strokes_topology, cam_param, 'lf')
% objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology, cam_param, 'scaleXYZ')
% % rotate3DAndSaveSVGFrames(strokes_topology, cam_param, '360');
% 
% [ strokes_topology, intersections, cam_param] = ...
% readReconstructionJsonFelix( fullfile(folder, 'student8_house_bestScore_full.json' ));
% [strokes_topology, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology, intersections, cam_param, []);
% folder_save = fullfile(foler_save_base, 'student8_house');
% %objectLFEffectAndSaveSVGFrames(strokes_topology, cam_param, 'lf')
% objectChangeScaleOneAxisAndSaveSVGFrames(strokes_topology, cam_param, 'scaleXYZ')
% % rotate3DAndSaveSVGFrames(strokes_topology, cam_param, '360');