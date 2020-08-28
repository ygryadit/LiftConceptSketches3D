function [folder_black, folder_color] = setupFOlder(folder_name)
    global folder_save;
    global save_color;


    folder_color = '';
    
    folder_svg = fullfile(folder_save, 'animation', folder_name);

    if ~exist(folder_svg, 'dir')
        mkdir(folder_svg);
    end

    folder_black = fullfile(folder_svg, 'black');
    if ~exist(folder_black , 'dir')
        mkdir(folder_black );
    end

    if save_color
        folder_color = fullfile(folder_svg, 'color');
        if ~exist(folder_color, 'dir')
            mkdir(folder_color);
        end
    end
end