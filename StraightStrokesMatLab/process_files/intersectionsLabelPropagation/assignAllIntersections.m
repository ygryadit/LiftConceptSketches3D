function [strokes_topology, intersections] = assignAllIntersections(folder_files, ...
                                                designer,...
                                                object,...
                                                view)

%% Data
% folder_files = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\_400';
% designer = 'student8';
% object = 'house';
% view = 'view1';
filename = [designer, '_', object, '_bestScore_full.json'];

%% Read file
fullfile(folder_files, designer, object, view, filename)
[ strokes_topology, intersections, ~] = ...
            readReconstructionJson( fullfile(folder_files, designer, object, view, filename));   
        

num_assigned = 0;
for i = 1:length(intersections)
   
   if  ~isempty(intersections(i).is_active)
        continue;
   end
   
   strokes_indices = intersections(i).strokes_indices;
   
   if (strokes_topology(strokes_indices(1)).primitive_type~=0) | ...
       (strokes_topology(strokes_indices(1)).primitive_type~=0)
       continue;
   end
   
   if (strokes_topology(strokes_indices(1)).depth_assigned==0) | ...
       (strokes_topology(strokes_indices(2)).depth_assigned==0)
       continue;
   end
   
   ic3d =             intersections(i).coordinates3D;
   
   stroke1points = readPointsData(strokes_topology(strokes_indices(1)));
   stroke2points = readPointsData(strokes_topology(strokes_indices(2)));
   
   length3D1 = strokes_topology(strokes_indices(1)).length3D;
   length3D2 = strokes_topology(strokes_indices(2)).length3D;
   
   belongs1 = checkInfIntersectionBelongToTheLine(stroke1points, length3D1,ic3d);
   belongs2 = checkInfIntersectionBelongToTheLine(stroke2points, length3D2,ic3d);
   
   intersections(i).is_active = belongs1 & belongs2;   
   num_assigned = num_assigned +1;
end

sketch_height = 512;
folder_save = fullfile(folder_files, designer, object, view);
save_as_svg_intersections_activation(strokes_topology,intersections, folder_save, sketch_height);
save_as_svg_intersections_activation_only_active(strokes_topology,intersections, folder_save, sketch_height);
save_as_svg_intersections_activation_only_non_active(strokes_topology,intersections, folder_save, sketch_height);

saveDrawingAsOBJ(strokes_topology, intersections, folder_save, 'propagatedIntersections');

end


function strokePoints = readPointsData(stroke)
    if ~isempty(stroke.points3D_clean)
        strokePoints = stroke.points3D_clean;
    else
        strokePoints = stroke.points3D;
    end
end

function belongs = checkInfIntersectionBelongToTheLine(points3D, length3D, intesetcion_coord3D)
    belongs = false;
    
    for i = 2:size(points3D,1)
        try
            distances_eq = pointLineDistance3D(intesetcion_coord3D, ...
                                               points3D(i,:),...
                                               points3D(i-1,:));
        catch e
           rethrow(e)
        end
         if distances_eq < 0.05*length3D                          
             belongs = true;
             return;
         end
    end
end

