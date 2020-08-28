function [folder_save, ...
 folder_save_imgs, ...
 filepath_sketch_json,...
 filepath_sketch_img] = setupFolderPaths(folder_save)

global datatset_name;        

if strcmp(datatset_name, 'OpenSketch')
    global folder_designer;
    global designer;
    global object_name;
    global view;

    folder_save = fullfile(folder_save, designer, object_name, view);

    if ~exist(folder_save, 'dir')
        mkdir(folder_save);
    end

    folder_save_imgs = fullfile(folder_save, 'images');
    if ~exist(folder_save_imgs, 'dir')
        mkdir(folder_save_imgs);
    end

    filepath_sketch_json = getFilepathSketchJson(   folder_designer, ...
                                                    designer,...
                                                    object_name,...
                                                    view,...                                            
                                                    'concept');
                                                
    filepath_sketch_img = fullfile(folder_designer,...
                                   designer,...
                                   object_name,...
                                   [view '_' 'concept' '_opaque.png']);



end

if strcmp(datatset_name, 'AnalyticDrawing')
    global folder_designer;
    global object_name;
    disp(folder_designer)


    folder_save = folder_designer
    folder_save_imgs = fullfile(folder_save, 'images');
    if ~exist(folder_save_imgs, 'dir')
        mkdir(folder_save_imgs);
    end
    filepath_sketch_json = fullfile(folder_designer,...
                                    [object_name '_clean.json'])
    filepath_sketch_img = fullfile(folder_designer,...
                                   [object_name '_clean.png'])
end

end