function [num_processed_strokes_skecthes,num_processed_ints_skecthes, T] = timingAnalysis()

folder1 = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\_400';

desingers = get_list_dirs_folder(folder1);
view = 'view1';

num_consistent_all = 0;
num_total = 0;

filenames = {};

times1 = [];
num_processed_strokes_skecthes = [];
num_processed_ints_skecthes =[];

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

        
        if ~exist(fullfile(folder1,designer, object, view, filename1), 'file') | ...
           ~exist(fullfile(folder1,designer, object, view, 'preformance.mat'), 'file')
            continue;
        end
        
        filenames{num_files} = filename1;
      
        
        [ strokes_topology_1, intersections_1, ~] = ...
            readReconstructionJson( fullfile(folder1,designer, object, view, filename1));   
        
        
        
        load(fullfile(folder1,designer, object, view, 'preformance.mat'), 'ellapsed_time');
        times1(num_files) = ellapsed_time/60;
        
        
        num_processed_strokes_skecthes(num_files) = sum([strokes_topology_1(:).primitive_type] == 0);
        num_processed_ints_skecthes(num_files) = sum([intersections_1(:).likely]);
        
        
        num_files = num_files +1;
    end
end


T = table(filenames', times1', num_processed_ints_skecthes', num_processed_strokes_skecthes');

[R1,p1] = corr([times1',num_processed_strokes_skecthes'],'Type','Spearman');
[R2,p2] = corr([times1',num_processed_ints_skecthes'],'Type','Spearman');

min(times1)
max(times1)
median(times1)
mean(times1)


end

function list = get_list_dirs_folder(folder)
desingers = dir(folder);
desingers = {desingers.name};
list = desingers(3:end);
end

function isIntersectionBetweenLines(intersection, strokes_topology_1)


end