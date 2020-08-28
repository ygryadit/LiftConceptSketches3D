function extractSVGS()

folder = 'C:\Users\yulia\Research\Data\sketches_json_first_viewpoint';

designers = dir(folder);
designers = {designers.name};
designers = designers(3:end);

for d = 1:length(designers )
    folder_d = fullfile(folder, designers{d});
    objects = dir(folder_d);
    objects = {objects.name};
    objects = objects(3:end);
    
    for o = 1:length(objects)
        launchOpenSketch(folder, folder_d, objects{o}, 'C:\Users\yulia\Research\DesignSketch3D\results\svgs')
    end
end

end