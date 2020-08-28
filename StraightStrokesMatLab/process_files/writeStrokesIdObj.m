folder_files = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\_400';
folder_curves = 'C:\Users\yulia\Research\DesignSketch3D\results\backproject_reconstruction_json\curves17052020\results\results_medium';
folder_save_obj = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\full_objs';
folder_save_json = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\full_json';
folder_save_svg = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\full_svg';

if ~exist(folder_save_json, 'dir')
    mkdir(folder_save_json);
end

if ~exist(folder_save_obj, 'dir')
    mkdir(folder_save_obj);
end

if ~exist(folder_save_svg, 'dir')
    mkdir(folder_save_svg);
end

global folder_save;
folder_save = folder_save_obj;

global sketch_height;
sketch_height = 512;

global designer;
global object;
global object_name;

files = dir(folder_curves);
files = files(~[files(:).isdir]);
files = {files.name};

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

   try
       [ strokes_topology_curves, intersections, cam_param] = ...
            readReconstructionJsonFelix( fullfile(folder_curves, filename));   
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
        strokes_topology_curves(ii).points2D = strokes_topology_rough(ii).points2D;
        strokes_topology_curves(ii).points2DOriginal = strokes_topology_rough(ii).points2DOriginal;
        strokes_topology_curves(ii).created = strokes_topology_rough(ii).created;
        strokes_topology_curves(ii).assigned = strokes_topology_rough(ii).assigned;
        strokes_topology_curves(ii).mean_pressure = strokes_topology_rough(ii).mean_pressure;
        strokes_topology_curves(ii).score = strokes_topology_rough(ii).score;
        strokes_topology_curves(ii).confidence = strokes_topology_rough(ii).confidence;
        
        if ~isnan(strokes_topology_curves(ii).points3D)
            strokes_topology_curves(ii).primitive_type = 0;
            strokes_topology_curves(ii).depth_assigned = true;
        end
    end
    
          
   [strokes_topology_curves, intersections, cam_param] = scaleAndCenterTheObject(strokes_topology_curves, intersections, cam_param, []);
   
   object_name = object;
   saveJSONReconstruction(strokes_topology_curves, intersections, cam_param, folder_save_json, 'rough');    
   saveDrawingAsOBJStrokesASObjects(strokes_topology_curves, folder_save_obj, 'rough');
    
   filename_svg = strrep(filename, '.json', '.svg');
   save_as_svg_sketch_original_rough(strokes_topology_curves, folder_save_svg, cam_param, filename_svg)
end


