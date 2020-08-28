function createAnimations(path, view)

if ~exist('path', 'var')
    path = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v_37_b\sketches';
end
if ~exist('view', 'var')
    view = 'view1';
end

desingers = dir(path);

desingers_names = {desingers.name};
desingers_names = desingers_names(3:end);
%desingers
for i = 1:length(desingers_names)
    designer = desingers_names{i};
    
    path_designer = fullfile(path, designer);
    
    objetcs = dir(path_designer);
    objetcs_names = {objetcs.name};
    objetcs_names = objetcs_names(3:end);
    
    %objects
    for o = 1:length(objetcs_names)
        object = objetcs_names{o};
        path_object = fullfile(path_designer, object, view);
        global folder_save;
        folder_save = path_object;
        
        folder_ = fullfile(folder_save, 'animation');
        if exist(folder_, 'dir')
            try
                rmdir(folder_, 's');
            catch

            end
        end
       try
        name = 'bestScore_full';
        createAnimationOneFile(path_object, designer, object, name);
%         
        name = 'highScore_full';
        createAnimationOneFile(path_object, designer, object, name);
        
        path_object = fullfile(path_designer, object, view, 'centered_reconstruction');
        name = 'final_full';
        createAnimationOneFile(path_object, designer, object, name);
       catch
       end
%         path_object = fullfile(path_designer, object, view, 'best_solutions');
%         
%         for i = 1:10
%             name = sprintf('highScore%04d_full', i);
%             createAnimationOneFile(path_object, designer, object, name);
%         end

    end

end





%confident 
%bestScore
%highScore

%10 best solutions:



end




