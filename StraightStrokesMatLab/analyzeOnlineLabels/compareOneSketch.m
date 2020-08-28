function [avg_val, interval_vals,...
         agreement, solution_labeling, solution_our, solution_designer,...
         percentage_agrees, ...
         avg_nd_agreemnt_match_designer, ...
         avg_nd_agreemnt_nmatch_designer,...
         agreement_designer,...
         mean_p] = compareOneSketch(folder, obj, folder_results, task)


%% List of all participants
participants = dir(folder);
participants = participants(3:end);
participants = participants(cat(1,participants.isdir));
participants = {participants.name}; 

%% Read all the labeling data
filename = 'intersections_active_0.json';
num_participants = length(participants)

labelings = [];
j = 1;
for i = 1:num_participants
   f = fullfile(folder, participants{i}, filename);
   try
        data = loadLabeling(f);
   catch
       continue;
   end
   labeling = struct2array(data.do_intersect);
   num_intrsctns = length(labeling);
   if i == 1
       labelings = zeros(num_participants,num_intrsctns);
   end
   if length(labeling) == size(labelings,2)
        labelings(j,:) = labeling;
        j = j +1;

   end
end

num_of_ones  = sum(labelings,1);
num_of_zeros = sum(~labelings,1);

[~,inds] = max([num_of_ones; num_of_zeros], [], 1);
agrement = [num_of_ones; num_of_zeros]./num_participants;
% linds = sub2ind(size(agrement),inds,1:length(inds));

agrement = max(agrement, [], 1);

solution = inds;
solution(solution == 2) = 0;
solution = logical(solution);

for i = 1:size(labelings)
    p(i) = sum(labelings(i,:) == solution)/length(solution);    
end
mean_p = mean(p);


%Compare to a designer:
[~, ...
 ~, ...
 ~, ...
 solution_designer,...
 agreement_designer] = compareSolution2Designer(solution, agrement, task);


% Load image:
img = imread(fullfile(folder, 'view1_concept_opaque.png'));
sketch      = readSketchJson(fullfile(folder, 'view1_concept.json'));
% sketch      = keepOnlyNonRemovedStrokes(sketch);


% Load intersections:
intersections = loadLabeling(fullfile(folder, 'intersections.json'));
[sketch, img, scale, intersections] = scaleCenterSketchIntersections(sketch, intersections, img);

close all;
fclose all;
f1 = figure(1);
imshow(img);
hold on;
intersections = selectSubsetIntersections(intersections);
folder_save = fullfile('E:\Projects\WiresProject\IntersectionsLabeling\results', obj);
if ~exist(folder_save, 'dir')
    mkdir(folder_save);
end
sketch_height = size(img,1);

plot(intersections.coordinates2D(solution,1), intersections.coordinates2D(solution,2), 'g*') 
plot(intersections.coordinates2D(~solution,1), intersections.coordinates2D(~solution,2), 'ro') 
legend('intersect', 'do not intersect');
saveas(f1, fullfile(folder, 'solution.png'));
text(intersections.coordinates2D(1:199,1)+2.5, intersections.coordinates2D(1:199,2)+2.5, string(round(agrement*100)/100.0))

colormap parula
colorbar
lim = caxis;
caxis([min(agrement) max(agrement)])



saveas(f1, fullfile(folder_save, 'solution_agrement.pdf'));
disp(folder)
interval_vals(1,1) =  mean(agrement(1:50));
interval_vals(1,2) =  mean(agrement(51:100));
interval_vals(1,3) =  mean(agrement(101:150));
interval_vals(1,4) =  mean(agrement(150:end));

avg_val = mean(agrement);

% disp(folder)
% disp(mean(agrement(1:50)))
% disp(mean(agrement(51:100)))
% disp(mean(agrement(101:150)))
% disp(mean(agrement(150:end)))
% 
% disp(mean(agrement()))

agrement_= agrement;
agrement_ = agrement_ - min(agrement_);
agrement_ = agrement_ / max(agrement_);



save_as_svg_straigt_strokes_intersect(sketch.strokes,intersections,solution, agrement_,folder_save,sketch_height,'intersections_active.svg')
save_as_svg_straigt_strokes_intersect(sketch.strokes,intersections,~solution, agrement_,folder_save,sketch_height,'intersections_nonactive.svg')



filepath = fullfile(folder_results, [obj '_bestScore_full.json'])

global designer;
global object_name;
   
   members = split(obj, '_');
   designer = members{1};
   
   object_parts = members(2:end);
   object_name = object_parts{1};
   for j = 2:length(object_parts)
      object_name = sprintf('%s_%s', object_name, object_parts{j});
   end
   
   
   folder_files = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\_400';
   
   
   
if ~exist(filepath, 'file')
    return;
end

[stroke_topology, intersections_rec] = assignAllIntersections(folder_files, ...
                                                designer,...
                                                object_name,...
                                                'view1');



% [stroke_topology, intersections_rec, cam_param] = readReconstructionJson(filepath);


%% Find intersections in the new data structure
% [sketch, img, scale, intersections] = scaleCenterSketchIntersections(sketch, intersections, img);
 map = mapStrokes(stroke_topology, sketch.strokes);
 
 intersections.strokes_indices = map(intersections.strokes_indices);
 map2 = toMerged(stroke_topology);
 for i = 1:size(intersections.strokes_indices,1)
    if ~isnan(intersections.strokes_indices(i,1))
        intersections.strokes_indices(i,1) = map2(intersections.strokes_indices(i,1));
    end
    if ~isnan(intersections.strokes_indices(i,2))
        intersections.strokes_indices(i,2) = map2(intersections.strokes_indices(i,2));
    end
 end
 
 inds_new = matchIntersections(intersections_rec, intersections, stroke_topology);

 
 ours_active = intersections.strokes_indices(:,1) == intersections.strokes_indices(:,2);
  
 inds = find(ours_active);
 intersections_coordinates2D_act = intersections.coordinates2D(inds,:);
 
 [agreement, solution_labeling, solution_our] = compareSolutions(solution, agrement, inds_new, intersections_rec, ours_active, intersections_coordinates2D_act, folder_save, stroke_topology);
 [agreement_designer, solution_designer, solution_our] = compareSolutions(solution_designer, agreement_designer, inds_new, intersections_rec, ours_active, intersections_coordinates2D_act, folder_save, stroke_topology);
 
    for i = 1:size(labelings,1)
        [~, labelings_(i,:), ~] = compareSolutions(labelings(i,:), agrement, inds_new, intersections_rec, ours_active, intersections_coordinates2D_act, folder_save, stroke_topology);
        p(i) = sum(labelings_(i,:) == solution_labeling)/length(solution_labeling);    
    end
    mean_p = mean(p);



   agrement_D_ND = (solution_designer == solution_labeling); %designer -- non-designer agrement
   percentage_agrees  = sum(agrement_D_ND)/length(solution_labeling);
   fprintf('percentage_agrees = %.2f\n', percentage_agrees);
   
   avg_nd_agreemnt_match_designer = mean(agreement(agrement_D_ND));
   avg_nd_agreemnt_nmatch_designer = mean(agreement(~agrement_D_ND));
end

function intersections_new = selectSubsetIntersections(intersections)

    num_intersections = size(intersections.coordinates2D,1);
    num_batch = 50;
    upper_num_inter = num_batch*4;
    if (num_intersections < upper_num_inter)
        intersections_new = intersections;
        return;
    end

    spacing = floor((num_intersections - upper_num_inter)/3.0);
    intersections_new.coordinates2D = zeros(upper_num_inter,2);
    intersections_new.strokes_indices = zeros(upper_num_inter,2);

    for j = 0:3
        for i = 0:(num_batch-1)            
            intersections_new.coordinates2D(i + j*num_batch+1,:)    = intersections.coordinates2D(i+j*num_batch+j*spacing+1,:);            
            intersections_new.strokes_indices(i + j*num_batch+1,:)  = intersections.strokes_indices(i+j*num_batch+j*spacing+1,:);
        end           
    end    
    
end