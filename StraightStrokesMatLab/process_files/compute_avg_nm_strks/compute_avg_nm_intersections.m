function T = compute_avg_nm_intersections(path_folder)
    num_intrsctns_vec = zeros(500, 1);
    num_lkl_intrsctns_vec = zeros(500, 1);
    
    num_sketches = 0;
    
    if ~exist('path_folder', 'var')
        path_folder = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\_400';
    end

    designers = getListOfFolders(path_folder);

    for d = 1:length(designers)
        path_designer = fullfile(path_folder, designers{d});
        objects = getListOfFolders(path_designer);
        for o = 1:length(objects)
            
            
            path_object = fullfile(path_designer, objects{o}, 'view1');
            path_stats = fullfile(path_object, 'intersections_stat.mat');
            if exist(path_stats, 'file')
                num_sketches = num_sketches+1;
                load(path_stats, 'num_intersections', 'num_likely_intersections');
                num_intrsctns_vec(num_sketches) = num_intersections;
                num_lkl_intrsctns_vec(num_sketches) = num_likely_intersections;
            end
        end
    end
    
    % Keep the data only:
    num_intrsctns_vec = num_intrsctns_vec(1:num_sketches);
    num_lkl_intrsctns_vec = num_lkl_intrsctns_vec(1:num_sketches);
    
    %Compute stats:
    avg_intrsctns = mean(num_intrsctns_vec);
    min_intrsctns = min(num_intrsctns_vec);
    max_intrsctns = max(num_intrsctns_vec);
    std_intrsctns = std(num_intrsctns_vec);
    
    avg_intrsctns_lkl = mean(num_lkl_intrsctns_vec);
    min_intrsctns_lkl = min(num_lkl_intrsctns_vec);
    max_intrsctns_lkl = max(num_lkl_intrsctns_vec);
    std_intrsctns_lkl = std(num_lkl_intrsctns_vec);
    
    T = table(avg_intrsctns,...
              min_intrsctns,...
              max_intrsctns,...
              std_intrsctns,...
              avg_intrsctns_lkl,...
              min_intrsctns_lkl,...
              max_intrsctns_lkl,...
              std_intrsctns_lkl);
end


function folders_names = getListOfFolders(path_folder)

folders_names = dir(path_folder);
folders_names = {folders_names.name};
folders_names = folders_names(3:end);
end