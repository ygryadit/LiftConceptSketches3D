% Statisctics for the section 7.4

function computePercentageOfDelayedStrokesWeakReduction()

%     folder_images = 'C:\Users\yulia\Research\Data\sketches_json_first_viewpoint';
    folder_files = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\_400';

%     global folder_save;
%     global sketch_height;
%     sketch_height = 512;

    designers = get_list_dirs_folder(folder_files);

    num_before = [];
    num_after = [];
    percentage_assigned_immediately_sketch = [];
    percentage_weak_affect_sketch = [];
    
    
    
    for d = 1:length(designers)
       close all
       designer = designers{d};
       objects = get_list_dirs_folder(fullfile(folder_files, designer));

       for o = 1:length(objects)
            object = objects{o};
            % Read sketches    
            filename1 = [designer '_' object '_bestScore_full.json'];
            
            if exist(fullfile(folder_files, designer, object, 'view1', filename1), 'file')
                [ strokes_topology, ~, ~] = ...
                    readReconstructionJson( fullfile(folder_files, designer, object, 'view1', filename1));   
            else
                continue;
            end
            [  inds_processed,...
            inds_assigned_immediately,...
            num_before_trim,...
            num_after_trim] = ...
                findIndicesProcessedStrokes(strokes_topology);
            
            [percentage_assigned_immediately_sketch(end+1),...
             percentage_weak_affect_sketch(end+1)] =...
               sketchStats( inds_processed,...
                            inds_assigned_immediately,...
                            num_before_trim,...
                            num_after_trim);
            num_before = [num_before, num_before_trim];
            num_after = [num_after, num_after_trim];
       end
        
    end

    mean_num_before = mean(num_before);
    mean_num_after = mean(num_after);
    mean_percentage_assigned_immediately_sketch = mean(percentage_assigned_immediately_sketch);
    mean_percentage_weak_affect_sketch = mean(percentage_weak_affect_sketch);

    
end


function [percentage_assigned_immediately, percentage_weak_affect] =...
               sketchStats( inds_processed,...
                            inds_assigned_immediately,...
                            num_before_trim,...
                            num_after_trim)
        
    percentage_assigned_immediately = length(inds_assigned_immediately)/...
        length(inds_processed);

    percentage_weak_affect= sum(num_before_trim ~= num_after_trim)/length(num_after_trim);

end
function [  inds_processed,...
            inds_assigned_immediately,...
            num_before_trim,...
            num_after_trim] = ...
                findIndicesProcessedStrokes(strokes_topology)
    
    inds_processed = [];
    inds_assigned_immediately = [];
    num_before_trim = [];
    num_after_trim = [];
    for i = 1:length(strokes_topology)
        if ~isempty(strokes_topology(i).created)
            inds_processed(end+1) = i;
            
            if (strokes_topology(i).created == strokes_topology(i).assigned)
                inds_assigned_immediately(end+1) = i;
            else
                num_before_trim(end+1) = strokes_topology(i).num_candidate_lines_before_trim;
                num_after_trim(end+1) = strokes_topology(i).num_candidate_lines_after_trim;
            end
        end
    end
end

function list = get_list_dirs_folder(folder)
desingers = dir(folder);
desingers = {desingers.name};
list = desingers(3:end);
end

function [object, designer, view] = parseObjectDesigner(filename)

filename_ = strrep(filename, '.json', '');
members = split(filename_,'_');
designer = members{1};

object_parts = members(2:(end-2));
object = object_parts{1};
for j = 2:length(object_parts)
  object = sprintf('%s_%s', object, object_parts{j});
end

view = 'view1';
       
end