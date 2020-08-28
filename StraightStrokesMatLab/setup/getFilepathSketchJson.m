% function getFilepathSketchJson(folder_designer, designer, object_name, view)
% 
% Returns the filepath of the json sketch representation.
% Input:
%   folder_designer -- path to the folder with subfolders of designers
%   designer        -- the name of designer, e.g. 'student1',
%                       'Professional1'.
%   object_name     -- the name of the object, e.g. 'vacuum cleaner' 
%   view            -- 'view1' or 'view2'
%   drawing_type    -- 'concept' or 'presentation', default 'concept'
% 
% Output:
%   filepath    -- filepath to json file with a sketch.
function filepath = getFilepathSketchJson(folder_designer, designer, object_name, view, drawing_type)
    if ~exist('drawing_type', 'var')
        drawing_type = 'concept';
    end
    
    filepath = fullfile(folder_designer,...
                        designer,...
                        object_name, ...
                        [view '_' drawing_type '.json']);
end