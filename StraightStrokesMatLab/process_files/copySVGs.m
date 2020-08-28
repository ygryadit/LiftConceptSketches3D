function copySVGs(folder_in, folder_out, view)
if ~exist('view', 'var')
    view = 'view1';
end


desingers = dir(folder_in);

desingers_names = {desingers.name};
desingers_names = desingers_names(3:end);
ni = 0;
nj = 0;
%desingers
for i = 1:length(desingers_names)
    designer = desingers_names{i};
    
    path_designer = fullfile(folder_in, designer);
    
    objetcs = dir(path_designer);
    objetcs_names = {objetcs.name};
    objetcs_names = objetcs_names(3:end);
    
    %objects
    for o = 1:length(objetcs_names)
        object = objetcs_names{o};
        path_svgs_in = fullfile(path_designer, object, view, 'lines_separation');
        path_svgs_out = fullfile(folder_out, designer, object, view, 'lines_separation');
        nj = nj+1;
        
        if ~exist(path_svgs_out, 'dir')
            mkdir(path_svgs_out)
        end
        
        file_straight_in = fullfile(path_svgs_in, 'lines.svg');
        file_straight_out = fullfile(path_svgs_out, 'lines.svg');
        if exist(file_straight_in, 'file')
            copyfile(file_straight_in, file_straight_out);
        end
        
        file_curves_in = fullfile(path_svgs_in, 'curves.svg');
        file_curves_out = fullfile(path_svgs_out, 'curves.svg');
        if exist(file_curves_in, 'file')
            copyfile(file_curves_in, file_curves_out);
            ni = ni + 1;
        end
        
    end
end

fprintf('Number files = %d out of %d\n', ni, nj);

end
