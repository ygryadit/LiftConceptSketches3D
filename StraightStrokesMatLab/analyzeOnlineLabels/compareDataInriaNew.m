function compareDataInriaNew()

folder_top = 'C:\Users\yulia\Research\DesignSketch3D\IntersectionsLabeling\data';
folder_results = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\reconstructions_json';

folders_tasks = dir(folder_top);
folders_tasks = folders_tasks(3:end);
folders_tasks = folders_tasks(cat(1,folders_tasks.isdir));
folders_tasks = {folders_tasks.name}; 
avg_val = zeros(3,1);
interval_vals = zeros(3,4);

for i = 1:3
    folder = fullfile(folder_top, folders_tasks{i});
    [avg_val(i), interval_vals(i,:),...
    agreement{i}, ...
    solution_labeling{i}, ...
    solution_our{i}, ...
    solution_designer{i},...
    percentage_agrees(i), ...
    avg_nd_agreemnt_match_designer(i), ...
    avg_nd_agreemnt_nmatch_designer(i),...
    agreement_designer{i},...
    mean_p] = compareOneSketch(folder,folders_tasks{i}, folder_results, folders_tasks{i});
end



interval_vals = interval_vals';

figure(1);
mx=max(interval_vals(:)) ;
colors = parula(3);
% plot([0.5,4.5], [avg_val(1) avg_val(1)], 'color', colors(1,:),'LineWidth',3)
% plot([0.5,4.5], [avg_val(2) avg_val(2)], 'color', colors(2,:),'LineWidth',3)
% plot([0.5,4.5], [avg_val(3) avg_val(3)], 'color', colors(3,:),'LineWidth',3)

hbar = bar(interval_vals,'FaceColor','flat');

for k = 1:size(interval_vals,2)
   hbar(k).CData = k;
end
% Set the axes YLIM (increaed wrt the max data value to have room for the
% label
ylim([0 mx*1.2])

% Get the XDATA
XDATA=get(hbar(1),'XData')';
% Define the vertical offset of the labels
ygap=mx*0.1;
% Loop over the bar's group
for i=1:length(hbar)
   % Get the YDATA of the i-th bar in each group
   YDATA=get(hbar(i),'YData')';
   % Loop over the groups
   for j=1:length(XDATA)
      % Get the XPOS of the j-th bar 
      xpos=XDATA(j);
      % Get the height of the bar and increment it with the offset
      ypos=YDATA(j)+ygap;
      % Define the labels
      labels=[num2str(YDATA(j),3)];
      % Add the labels
      t = text(xpos+hbar(i).XOffset,ypos,labels,'Color','k','HorizontalAlignment','center','Rotation',90)
   end
end 

legend(folders_tasks,'location','southoutside')


%%
agreement_all = cat(2,agreement{:});
agreement_designer_all = cat(2,agreement_designer{:});

solution_labeling_all = cat(2,solution_labeling{:});
solution_our_all = cat(1,solution_our{:})';
solution_designer_all = cat(2,solution_designer{:});

 
 mask_agreement_75 =  agreement_all > 0.75;
 mask_agreement_50 =  (agreement_all > 0.5) & agreement_all < 0.75;
 
 percentage_same = sum((solution_our_all == solution_labeling_all))./length(solution_labeling_all);
 
 percentage_same_d = sum((solution_our_all == solution_designer_all))./length(solution_designer_all);
 
 
 percentage_same_75 = sum((solution_our_all(mask_agreement_75) == ...
                        solution_labeling_all(mask_agreement_75)))./...
                        length(solution_labeling_all(mask_agreement_75));
                    
 percentage_same_50 = sum((solution_our_all(mask_agreement_50) == ...
                        solution_labeling_all(mask_agreement_50)))./...
                        length(solution_labeling_all(mask_agreement_50));    
                    
                    
 fprintf('Full = %.3f, 75 = %.3f, 50 = %.3f\n', ...
            percentage_same, ...
            percentage_same_75, ...
            percentage_same_50)       
        
fprintf('Full designer= %.3f \n', ...
            percentage_same_d)    
        
fprintf('mean agreement designers nds= %.3f \n', ...
            mean(percentage_agrees))            
        
fprintf('mean avg_nd_agreemnt_match_designer= %.3f \n', ...
            mean(avg_nd_agreemnt_match_designer))            
        
fprintf('mean avg_nd_agreemnt_nmatch_designer= %.3f \n', ...
            mean(avg_nd_agreemnt_nmatch_designer))            
end

