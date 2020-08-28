function createAnimationOneFile(path, designer, object, name)

  filepath = getFilePath(path, designer, object, name);
  
  if ~exist(filepath, 'file')
      return;
  end
  
  [ strokes_topology, ~, cam_param] = ...
            readReconstructionJson( filepath );
       
global filepath_sketch_json;
global filepath_sketch_img;
folder_designer = 'C:\Users\yulia\Research\Data\sketches_json_first_viewpoint\';

   filepath_sketch_json = getFilepathSketchJson(   folder_designer, ...
                                                    designer,...
                                                    object,...
                                                    'view1',...                                            
                                                    'concept');
                                                
    filepath_sketch_img = fullfile(folder_designer,...
                                   designer,...
                                   object,...
                                   ['view1' '_' 'concept' '_opaque.png']);
                               
  [~, ~] = loadSketch();

  rotate3DAndSaveSVGFrames(strokes_topology, cam_param, name);
end


function filepath = getFilePath(path, designer, object, name)
    
    filepath = fullfile(path, ...
                sprintf('%s_%s_%s.json', designer, object, name) ...
                );

end