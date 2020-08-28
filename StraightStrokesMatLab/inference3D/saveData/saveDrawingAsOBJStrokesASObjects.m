function saveDrawingAsOBJStrokesASObjects(strokes_topology, folder_save, fname)
global designer;
global object_name;

global SAVE_SVGS;
if ~SAVE_SVGS
    return;
end

file_save_obj = fullfile(folder_save, ...
                         sprintf('%s_%s_%s.obj', ...
                                    designer, ...
                                    object_name, ...
                                    fname))
fid = fopen(file_save_obj , 'w');

str = sprintf('g [%s_%s]\n', designer, object_name); 
fwrite(fid, str);  
    
vj = 0;
for i = 1:length(strokes_topology)    
    
    if ~strokes_topology(i).depth_assigned
        continue;
    end
    
    
    vj = printOBJLines(1, i, strokes_topology, fid, vj);
    
end



end

function vj = printOBJLines(num_strokes, indices, strokes_topology, fid, vj)
    global ZUP;
    
    str = sprintf('o [%d]\n', indices); 
    fwrite(fid, str);  
    
    
    for k = 1:num_strokes
       i = indices(k);
       vj = vj + 1;
       for j = 1:size(strokes_topology(i).points3D,1)
           str = sprintf('v %.3f %.3f %.3f\n', strokes_topology(i).points3D(j,1), -strokes_topology(i).points3D(j,3), strokes_topology(i).points3D(j,2)); 
           fwrite(fid, str);  
       end 
       for j = 1:(size(strokes_topology(i).points3D,1)-1)
           str = sprintf('l %d %d\n', vj, vj+1); 
           vj = vj + 1;
           fwrite(fid, str); 
       end
    end
end