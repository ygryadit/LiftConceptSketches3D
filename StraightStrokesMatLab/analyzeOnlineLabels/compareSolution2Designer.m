function [percentage_agrees, avg_nd_agreemnt_match_designer, avg_nd_agreemnt_nmatch_designer, solution, agrement] = compareSolution2Designer(solution_ref, agrement_ref, task)
    folder_designer = 'C:\Users\yulia\Research\DesignSketch3D\IntersectionsLabeling\data_designer';
    folder = fullfile(folder_designer, task);
    
    %% List of all designers
    participants = dir(folder);
    participants = participants(3:end);
    participants = participants(cat(1,participants.isdir));
    participants = {participants.name};
    
    %% Read all the labeling data
    filename = 'intersections_active_0.json';
%     participants = participants(2);
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
    sum(agrement==1)/length(agrement)
    mean(agrement)
    solution = inds;
    solution(solution == 2) = 0;
    solution = logical(solution);
    
   %% Evalaute designers agrement on a solution with average solution 
   agrement_D_ND = (solution == solution_ref); %designer -- non-designer agrement
   percentage_agrees  = sum(agrement_D_ND)/length(solution_ref);
   fprintf('percentage_agrees = %.2f\n', percentage_agrees);
   
   avg_nd_agreemnt_match_designer = mean(agrement_ref(agrement_D_ND));
   avg_nd_agreemnt_nmatch_designer = mean(agrement_ref(~agrement_D_ND));
   
end