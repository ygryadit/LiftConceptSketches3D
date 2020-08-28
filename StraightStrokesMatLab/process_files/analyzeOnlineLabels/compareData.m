function compareData()

folder_object = 'Z:\WiresProject\Data\grth_labeling\student8\house';

f1 = fullfile(folder_object, 'Adrien\intersections_active_0.json');
f2 = fullfile(folder_object, 'Chenxi\intersections_active_0.json');

data1 =  loadLabeling(f1);
data2 =  loadLabeling(f2);

fields1 = fieldnames(data1.do_intersect);
fields2 = fieldnames(data2.do_intersect);

labeling1 = struct2array(data1.do_intersect)
labeling2 = struct2array(data2.do_intersect)

num_differ = sum(labeling1 ~= labeling2)
percentage_different = num_differ/length(labeling1)*100

mask_differ = (labeling1 ~= labeling2);

% Load image:
img = imread(fullfile(folder_object, 'view1_concept_opaque.png'));

% Load intersections:
intersections = loadLabeling(fullfile(folder_object, 'intersections.json'));

%%
figure;
[ha, pos] = tight_subplot(1,3,[.01 .03],[.1 .01],[.01 .01]) 

axes(ha(1))
imshow(img);
hold on;

plot(intersections.coordinates2D(labeling1,1), intersections.coordinates2D(labeling1,2), 'g*') 
plot(intersections.coordinates2D(~labeling1,1), intersections.coordinates2D(~labeling1,2), 'r*') 

legend('active', 'non-active');
title('User 1')

%%
axes(ha(2))
imshow(img);
hold on;

plot(intersections.coordinates2D(labeling2,1), intersections.coordinates2D(labeling2,2), 'g*') 
plot(intersections.coordinates2D(~labeling2,1), intersections.coordinates2D(~labeling2,2), 'r*') 

legend('active', 'non-active');
title('User 2')

%%
axes(ha(3))
imshow(img);
hold on;


plot(intersections.coordinates2D(~mask_differ,1), intersections.coordinates2D(~mask_differ,2), 'g*') 
plot(intersections.coordinates2D(mask_differ,1), intersections.coordinates2D(mask_differ,2), 'r*') 

legend('identical', 'differ');
title(sprintf('Differencies: %.2f percent of intersections label differ', percentage_different));

end

function data = loadLabeling(filepath)
   
   fid = fopen(filepath);
   raw = fread(fid,inf);
   str = char(raw');
   fclose(fid);
   data = jsondecode(str);
   
end