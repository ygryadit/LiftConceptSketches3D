lists;

%% Read json data for original strokes:
% keep only starigth stroke and those that have imformation about the
% assignement.
path_json = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\full_json';
percentage_delayed_dave_14 = [];
percentage_delayed_jerry_14 = [];

percentage_delayed_good_dave_14 = [];
percentage_delayed_good_jerry_14 = [];


for i = 1:size(list_designers_paper_14,1)
    desinger{i} = list_designers_paper_14{i,1};
    object{i} = list_designers_paper_14{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough_full.json']
    
   [strokes_topology_rough, ~, ~] = ...
            readReconstructionJson( fullfile(path_json, filename)); 
   
    % dave
    
   set_removed = strokes_numbers_dave_removed{i};
    
   pcl = computePercentageDelayedm(set_removed, strokes_topology_rough);
   percentage_delayed_dave_14(i) = pcl;
   
   percentage_delayed_good_dave_14(i) = computePercentageDelayedm(strokes_numbers_dave{i}, strokes_topology_rough);
   
   
   set_removed = strokes_numbers_jerry_removed{i};
   
   pcl = computePercentageDelayedm(set_removed, strokes_topology_rough);
   percentage_delayed_jerry_14(i) = pcl;
   percentage_delayed_good_jerry_14(i) = computePercentageDelayedm(strokes_numbers_jerry{i}, strokes_topology_rough);
end

fprintf('Mean_percentage_delayed_dave = %.2f\n', mean(percentage_delayed_dave_14(~isnan(percentage_delayed_dave_14))));
fprintf('Median_percentage_delayed_dave = %.2f\n', median(percentage_delayed_dave_14(~isnan(percentage_delayed_dave_14))));

fprintf('Mean_percentage_delayed_jerry = %.2f\n', mean(percentage_delayed_jerry_14(~isnan(percentage_delayed_jerry_14))));
fprintf('Median_percentage_delayed_jerry = %.2f\n', median(percentage_delayed_jerry_14(~isnan(percentage_delayed_jerry_14))));

maks_average = find(~isnan(percentage_delayed_dave_14) & ~isnan(percentage_delayed_jerry_14));
maks_dave = find(~isnan(percentage_delayed_dave_14) & isnan(percentage_delayed_jerry_14));
maks_jerry = find(isnan(percentage_delayed_dave_14) & ~isnan(percentage_delayed_jerry_14));

percentage_delayed_14 = NaN*ones(size(percentage_delayed_dave_14));
percentage_delayed_14(maks_average) = 0.5*(percentage_delayed_dave_14(maks_average)+percentage_delayed_jerry_14(maks_average));
percentage_delayed_14(maks_dave) = percentage_delayed_dave_14(maks_dave);
percentage_delayed_14(maks_jerry) = percentage_delayed_jerry_14(maks_jerry);

fprintf('Mean_percentage_delayed_14 = %.2f\n', mean(percentage_delayed_14(~isnan(percentage_delayed_14))));
fprintf('Median_percentage_delayed_14 = %.2f\n', median(percentage_delayed_14(~isnan(percentage_delayed_14))));

percentage_delayed_good_14 = combineLabels(percentage_delayed_good_dave_14, percentage_delayed_good_jerry_14);

%% supp
percentage_delayed_yuan_supp = [];
percentage_delayed_jerry_supp = [];

percentage_delayed_good_jerry_supp = [];
percentage_delayed_good_yuan_supp = [];

for i = 1:size(list_designers_supplemental,1)
    desinger{i} = list_designers_supplemental{i,1};
    object{i} = list_designers_supplemental{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough_full.json'];
    
   [strokes_topology_rough, ~, ~] = ...
            readReconstructionJson( fullfile(path_json, filename)); 
   
    % yuan
    
   set_removed = strokes_numbers_yuan_removed_supp{i};
    
   pcl = computePercentageDelayedm(set_removed, strokes_topology_rough);
   percentage_delayed_yuan_supp(i) = pcl;
   percentage_delayed_yuan_dave_supp(i) = computePercentageDelayedm(strokes_numbers_yuan_supp{i}, strokes_topology_rough);
   % jerry
   set_removed = strokes_numbers_jerry_removed_supp{i};
   
   pcl = computePercentageDelayedm(set_removed, strokes_topology_rough);   
   percentage_delayed_jerry_supp(i) = pcl;
   
   percentage_delayed_good_jerry_supp(i) = computePercentageDelayedm(strokes_numbers_jerry_supp{i}, strokes_topology_rough);
end

fprintf('Mean_percentage_delayed_yuan_supp = %.2f\n', mean(percentage_delayed_yuan_supp(~isnan(percentage_delayed_yuan_supp))));
fprintf('Median_percentage_delayed_yuan_supp = %.2f\n', median(percentage_delayed_yuan_supp(~isnan(percentage_delayed_yuan_supp))));

fprintf('Mean_percentage_delayed_jerry_supp = %.2f\n', mean(percentage_delayed_jerry_supp(~isnan(percentage_delayed_jerry_supp))));
fprintf('Median_percentage_delayed_jerry_supp = %.2f\n', median(percentage_delayed_jerry_supp(~isnan(percentage_delayed_jerry_supp))));

maks_average = find(~isnan(percentage_delayed_yuan_supp) & ~isnan(percentage_delayed_jerry_supp));
maks_dave = find(~isnan(percentage_delayed_yuan_supp) & isnan(percentage_delayed_jerry_supp));
maks_jerry = find(isnan(percentage_delayed_yuan_supp) & ~isnan(percentage_delayed_jerry_supp));

percentage_delayed_supp = NaN*ones(size(percentage_delayed_yuan_supp));
percentage_delayed_supp(maks_average) = 0.5*(percentage_delayed_yuan_supp(maks_average)+percentage_delayed_jerry_supp(maks_average));
percentage_delayed_supp(maks_dave) = percentage_delayed_yuan_supp(maks_dave);
percentage_delayed_supp(maks_jerry) = percentage_delayed_jerry_supp(maks_jerry);

percentage_delayed_good_supp = combineLabels(percentage_delayed_good_jerry_supp, percentage_delayed_yuan_dave_supp);

v_supp = percentage_delayed_supp(~isnan(percentage_delayed_supp));
v14 = percentage_delayed_14(~isnan(percentage_delayed_14));

fprintf('Mean_percentage_delayed_supp = %.2f\n', mean(v_supp));
fprintf('Median_percentage_delayed_supp = %.2f\n', median(percentage_delayed_supp(~isnan(percentage_delayed_supp))));

fprintf('Mean_percentage_delayed_all = %.2f\n', mean( [v_supp  v14]));
fprintf('Median_percentage_delayed_all = %.2f\n', median([v_supp  v14]  ));

v_supp = percentage_delayed_good_supp(~isnan(percentage_delayed_good_supp));
v14 = percentage_delayed_good_14(~isnan(percentage_delayed_good_14));

fprintf('Mean_percentage_delayed_all = %.2f\n', mean( [v_supp  v14]));
fprintf('Median_percentage_delayed_all = %.2f\n', median([v_supp  v14]  ));


