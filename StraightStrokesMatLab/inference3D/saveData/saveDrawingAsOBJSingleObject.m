function saveDrawingAsOBJSingleObject(strokes_topology, folder_save, fname)
global designer;
global object_name;

global SAVE_SVGS;
if ~SAVE_SVGS
    return;
end

file_save_obj = fullfile(folder_save, ...
                         sprintf('%s_%s_%s_single_object.obj', ...
                         designer,...
                         object_name,...
                         fname));
                     
fid = fopen(file_save_obj , 'w');

str = sprintf('o [%s_%s]\n', designer, object_name); 
fwrite(fid, str);  
       
ind_consider = find(cat(1,strokes_topology.primitive_type)==0);

group_1 = find(cat(1,strokes_topology(ind_consider).line_group) == 1  & cat(1,strokes_topology(ind_consider).depth_assigned) == 1);
group_2 = find(cat(1,strokes_topology(ind_consider).line_group) == 2  & cat(1,strokes_topology(ind_consider).depth_assigned) == 1);
group_3 = find(cat(1,strokes_topology(ind_consider).line_group) == 3  & cat(1,strokes_topology(ind_consider).depth_assigned) == 1);
group_4 = find(cat(1,strokes_topology(ind_consider).line_group) == 4  & cat(1,strokes_topology(ind_consider).depth_assigned) == 1);
group_5 = find(cat(1,strokes_topology(ind_consider).line_group) == 5  & cat(1,strokes_topology(ind_consider).depth_assigned) == 1);


num_strokes_1 = length(group_1);
num_strokes_2 = length(group_2);
num_strokes_3 = length(group_3);
num_strokes_4 = length(group_4);
num_strokes_5 = length(group_5);

group_1 = ind_consider(group_1);
group_2 = ind_consider(group_2);
group_3 = ind_consider(group_3);
group_4 = ind_consider(group_4);
group_5 = ind_consider(group_5);

fwrite(fid, str);  
vj = 0;
vj = printOBJLines(num_strokes_1, group_1, strokes_topology, fid, vj);

vj = printOBJLines(num_strokes_2, group_2, strokes_topology, fid, vj);

vj = printOBJLines(num_strokes_3, group_3, strokes_topology, fid, vj);

vj = printOBJLines(num_strokes_4, group_4, strokes_topology, fid, vj);

vj = printOBJLines(num_strokes_5, group_5, strokes_topology, fid, vj);




end

function vj = printOBJLines(num_strokes, indices, strokes_topology, fid, vj)
    global ZUP
    if ~ZUP
        signz = -1;
    else
        signz = 1; 
    end
    s = 10;
    for k = 1:num_strokes
       i = indices(k);
       vj = vj + 1;
       for j = 1:size(strokes_topology(i).points3D,1)
           str = sprintf('v %.3f %.3f %.3f\n', strokes_topology(i).points3D(j,1)*s, ...
                                              -strokes_topology(i).points3D(j,3)*s*signz,...
                                              strokes_topology(i).points3D(j,2)*s); 
           fwrite(fid, str);  
       end 
       for j = 1:(size(strokes_topology(i).points3D,1)-1)
           str = sprintf('l %d %d\n', vj, vj+1); 
           vj = vj + 1;
           fwrite(fid, str); 
       end
    end
end