lists;

%% Read json data for original strokes:
% keep only starigth stroke and those that have imformation about the
% assignement.
path_json = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\full_json';
delayes_14 =[]

for i = 1:size(list_designers_paper_14,1)
    desinger{i} = list_designers_paper_14{i,1};
    object{i} = list_designers_paper_14{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough_full.json']
    
   [strokes_topology_rough, ~, ~] = ...
            readReconstructionJson( fullfile(path_json, filename)); 
   
    % dave
    
   set_removed = strokes_numbers_dave_removed{i};
    
   pcl = computeDelay(set_removed, strokes_topology_rough);
   if ~isnan(pcl)
        delayes_14 = [delayes_14; pcl];
   end
   
   set_removed = strokes_numbers_jerry_removed{i};
   
   pcl = computeDelay(set_removed, strokes_topology_rough);
   if ~isnan(pcl)
        delayes_14 = [delayes_14; pcl];
   end
end

% fprintf('Mean_percentage_assigned_last_dave = %.2f\n', mean(percentage_assigned_last_dave_14(~isnan(percentage_assigned_last_dave_14))));
% fprintf('Median_percentage_assigned_last_dave = %.2f\n', median(percentage_assigned_last_dave_14(~isnan(percentage_assigned_last_dave_14))));
% 
% fprintf('Mean_percentage_assigned_last_jerry = %.2f\n', mean(percentage_assigned_last_jerry_14(~isnan(percentage_assigned_last_jerry_14))));
% fprintf('Median_percentage_assigned_last_jerry = %.2f\n', median(percentage_assigned_last_jerry_14(~isnan(percentage_assigned_last_jerry_14))));
% 
% maks_average = find(~isnan(percentage_assigned_last_dave_14) & ~isnan(percentage_assigned_last_jerry_14));
% maks_dave = find(~isnan(percentage_assigned_last_dave_14) & isnan(percentage_assigned_last_jerry_14));
% maks_jerry = find(isnan(percentage_assigned_last_dave_14) & ~isnan(percentage_assigned_last_jerry_14));
% 
% percentage_assigned_last_14 = NaN*ones(size(percentage_assigned_last_dave_14));
% percentage_assigned_last_14(maks_average) = 0.5*(percentage_assigned_last_dave_14(maks_average)+percentage_assigned_last_jerry_14(maks_average));
% percentage_assigned_last_14(maks_dave) = percentage_assigned_last_dave_14(maks_dave);
% percentage_assigned_last_14(maks_jerry) = percentage_assigned_last_jerry_14(maks_jerry);

fprintf('Mean_percentage_assigned_last_14 = %.2f\n', mean(delayes_14));
fprintf('Median_percentage_assigned_last_14 = %.2f\n', median(delayes_14));


%% supp
percentage_assigned_last_yuan_supp = [];
percentage_assigned_last_jerry_supp = [];
delayes_supp =[];
for i = 1:size(list_designers_supplemental,1)
    desinger{i} = list_designers_supplemental{i,1};
    object{i} = list_designers_supplemental{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough_full.json'];
    
   [strokes_topology_rough, ~, ~] = ...
            readReconstructionJson( fullfile(path_json, filename)); 
   
    % yuan
    
   set_removed = strokes_numbers_yuan_removed_supp{i};
    
   pcl = computeDelay(set_removed, strokes_topology_rough);
   if ~isnan(pcl)
        delayes_supp = [delayes_supp; pcl];
   end
   
   
   % jerry
   set_removed = strokes_numbers_jerry_removed_supp{i};
   
   pcl = computeDelay(set_removed, strokes_topology_rough);
   if ~isnan(pcl)
        delayes_supp = [delayes_supp; pcl];
   end
end

% fprintf('Mean_percentage_assigned_last_yuan_supp = %.2f\n', mean(percentage_assigned_last_yuan_supp(~isnan(percentage_assigned_last_yuan_supp))));
% fprintf('Median_percentage_assigned_last_yuan_supp = %.2f\n', median(percentage_assigned_last_yuan_supp(~isnan(percentage_assigned_last_yuan_supp))));
% 
% fprintf('Mean_percentage_assigned_last_jerry_supp = %.2f\n', mean(percentage_assigned_last_jerry_supp(~isnan(percentage_assigned_last_jerry_supp))));
% fprintf('Median_percentage_assigned_last_jerry_supp = %.2f\n', median(percentage_assigned_last_jerry_supp(~isnan(percentage_assigned_last_jerry_supp))));
% 
% maks_average = find(~isnan(percentage_assigned_last_yuan_supp) & ~isnan(percentage_assigned_last_jerry_supp));
% maks_dave = find(~isnan(percentage_assigned_last_yuan_supp) & isnan(percentage_assigned_last_jerry_supp));
% maks_jerry = find(isnan(percentage_assigned_last_yuan_supp) & ~isnan(percentage_assigned_last_jerry_supp));
% 
% percentage_assigned_last_supp = NaN*ones(size(percentage_assigned_last_yuan_supp));
% percentage_assigned_last_supp(maks_average) = 0.5*(percentage_assigned_last_yuan_supp(maks_average)+percentage_assigned_last_jerry_supp(maks_average));
% percentage_assigned_last_supp(maks_dave) = percentage_assigned_last_yuan_supp(maks_dave);
% percentage_assigned_last_supp(maks_jerry) = percentage_assigned_last_jerry_supp(maks_jerry);
% 
% 
% v_supp = percentage_assigned_last_supp(~isnan(percentage_assigned_last_supp));
% v14 = percentage_assigned_last_14(~isnan(percentage_assigned_last_14));

fprintf('Mean_percentage_assigned_last_supp = %.2f\n', mean(delayes_supp));
fprintf('Median_percentage_assigned_last_supp = %.2f\n', median(delayes_supp));

fprintf('Mean_percentage_assigned_last_all = %.2f\n', mean( [delayes_supp;  delayes_14]));
fprintf('Median_percentage_assigned_last_all = %.2f\n', median([delayes_supp;  delayes_14]  ));
