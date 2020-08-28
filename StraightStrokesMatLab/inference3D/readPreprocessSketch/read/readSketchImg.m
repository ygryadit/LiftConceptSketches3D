function img = readSketchImg(filepath_sketch_img, READ)
   global SHOW_FIGS_PREPROCESS;
   global SHOW_FIGS;
   
   if ~exist('READ', 'var')
       READ = false;
   end
   
   if ~(READ || SHOW_FIGS_PREPROCESS || SHOW_FIGS )
       img = [];
       return;
   end
       
    if exist(filepath_sketch_img, 'file')
            img = imread(filepath_sketch_img);   
    else
        % Convert transparent version to opaque and save the opaque version        
        filename = strrep(filepath_sketch_img,  '_opaque.png', '.png');
        
        if exist(filename, 'file')
            [~,~,img] = imread(filename);   
            img = 255-img;

            if size(img,3) == 1
                img = repmat(img, 1,1,3);
            end
            imwrite(img, filepath_sketch_img);                       
        end
    end
    if size(img,3) == 1
            img = repmat(img, 1,1,3);
    end
end