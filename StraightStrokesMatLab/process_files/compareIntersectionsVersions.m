function [num_consistent_all,num_total, T] = compareIntersectionsVersions(folder1,folder2)

% folder1 = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\_400';
% % folder2 = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_all_inter_400\_400';
% folder2 = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400_no_merge\_400';


folder1 = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v51\_400';
folder2 = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v51\_400_now_weak_rejection';

desingers = get_list_dirs_folder(folder1);
view = 'view1';


num_consistent_all = 0;
num_total = 0;

filenames = {};
num_consistent_sketches = [];
times1 = [];
times2 = [];

num_files = 1;
% Go over designers
for d = 1:length(desingers)
    designer = desingers{d};
    objects = get_list_dirs_folder(fullfile(folder1, desingers{d}));
    % Go over objects
    for o = 1:length(objects)
        object = objects{o};
        % Read sketches    
        filename1 = [designer '_' object '_bestScore_full.json'];
        filename2 = [designer '_' object '_bestScore_full.json'];
        
        if ~exist(fullfile(folder1,designer, object, view, filename1), 'file') | ...
           ~exist(fullfile(folder2,designer, object, view, filename2), 'file') | ...
           ~exist(fullfile(folder1,designer, object, view, 'preformance.mat'), 'file') | ...
           ~exist(fullfile(folder2,designer, object, view, 'preformance.mat'), 'file')
            continue;
        end
        
        filenames{num_files} = filename1;
      
        
        [ strokes_topology_1, intersections_1, ~] = ...
            readReconstructionJson( fullfile(folder1,designer, object, view, filename1));   
        fullfile(folder2,designer, object, view, filename2)
        [ strokes_topology_2, intersections_2, ~] = ...
            readReconstructionJson( fullfile(folder2,designer, object, view, filename2));   
        
        
        load(fullfile(folder1,designer, object, view, 'preformance.mat'), 'ellapsed_time');
        times1(num_files) = ellapsed_time;
        load(fullfile(folder2,designer, object, view, 'preformance.mat'), 'ellapsed_time');
        times2(num_files) = ellapsed_time;
        
        
        num_consistent_sketch = 0;
        num_total_sketch = 0;
            
        % Go over the intersection is the first sketch
        for i = 1:length(intersections_1)
            if isempty(intersections_1(i).is_active)
                continue;
            end
            
            if intersections_1(i).is_active == intersections_2(i).is_active
                num_consistent_sketch = num_consistent_sketch + 1;
            end
            num_total_sketch = num_total_sketch + 1;
            % Compare labels
                    
        end
        num_consistent_all = num_consistent_all + num_consistent_sketch;
        num_total = num_total + num_total_sketch;
        
        num_consistent_sketches(num_files) = num_consistent_sketch/num_total_sketch;
        num_files = num_files + 1;
    end
end

T = table(filenames', num_consistent_sketches', times1', times2');
save('C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\no_merge.mat', 'num_consistent_all', 'num_total', 'T')
end

function list = get_list_dirs_folder(folder)
desingers = dir(folder);
desingers = {desingers.name};
list = desingers(3:end);
end

function isIntersectionBetweenLines(intersection, strokes_topology_1)


end