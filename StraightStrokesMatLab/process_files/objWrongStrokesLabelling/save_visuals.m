lists;

%% Read json data for original strokes:
% keep only starigth stroke and those that have imformation about the
% assignement.
path_json = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\full_json';


for i = 1:size(list_designers_paper_14,1)
    desinger{i} = list_designers_paper_14{i,1};
    object{i} = list_designers_paper_14{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough_full.json']
    
    [strokes_topology_rough, ~, cam_param] = ...
            readReconstructionJson( fullfile(path_json, filename)); 
   
 
   % dave 
   set_removed = strokes_numbers_dave_removed{i};
   R = rotz_mine(0);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_dave_view_1.svg'], cam_param, set_removed,R);
   R = rotz_mine(15);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_dave_view_2.svg'], cam_param, set_removed,R);
   R = rotz_mine(-15);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_dave_view_3.svg'], cam_param, set_removed,R);
   
   set_removed = strokes_numbers_jerry_removed{i};
   R = rotz_mine(0);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_jerry_view_1.svg'], cam_param, set_removed, R);
   R = rotz_mine(15);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_jerry_view_2.svg'], cam_param, set_removed, R);
   R = rotz_mine(-15);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_jerry_view_3.svg'], cam_param, set_removed,R);
end




%% supp
percentage_assigned_last_yuan_supp = [];
percentage_assigned_last_jerry_supp = [];

for i = 1:size(list_designers_supplemental,1)
    desinger{i} = list_designers_supplemental{i,1};
    object{i} = list_designers_supplemental{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough_full.json'];
    
   [strokes_topology_rough, ~, cam_param] = ...
            readReconstructionJson( fullfile(path_json, filename)); 
   
    % yuan
    
   set_removed = strokes_numbers_yuan_removed_supp{i};
   R = rotz_mine(0);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_yuan_view_1.svg'], cam_param, set_removed, R);
   R = rotz_mine(15);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_yuan_view_2.svg'], cam_param, set_removed, R);
   R = rotz_mine(-15);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_yuan_view_3.svg'], cam_param, set_removed, R);
   % jerry
   set_removed = strokes_numbers_jerry_removed_supp{i};
   R = rotz_mine(0);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_jerry_view_1.svg'], cam_param, set_removed, R);
   R = rotz_mine(15);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_jerry_view_2.svg'], cam_param, set_removed, R);
   R = rotz_mine(-15);
   save_as_svg_sketch_visual_bad_strokes(strokes_topology_rough, folder_save, [desinger{i} '_' object{i} '_jerry_view_3.svg'], cam_param, set_removed, R);
end



