clear all
list_designers_paper_14 = [ {'designer2'}, {'bookshelves'};
    {'designer2'}, {'cabinet_02'};
    {'designer2'}, {'chair_04'};
    {'designer2'}, {'printer_02'};
    {'Prof2task2'}, {'printer_02'};
    {'Professional4'}, {'tubes'};
    {'Professional6'}, {'house'};
    {'Professional6'}, {'mouse'};
    {'Professional6'}, {'shampoo_bottle'};
    {'Professional6'}, {'wobble_surface'};
    {'student1'}, {'vacuum_cleaner'};
    {'student7'}, {'wobble_surface'};
    {'student8'}, {'mouse'};
    {'student9'}, {'house'};
];

list_designers_supplemental = [ {'designer2'}, {'armchair_02'};
    {'designer2'}, {'cabinet_01'};
    {'designer2'}, {'guitar_01'};
    {'designer2'}, {'lamp'};
    {'Prof2task2'}, {'cabinet_01'};
    {'Prof2task2'}, {'guitar_01'};
    {'Professional1'}, {'house'};
    {'Professional1'}, {'tubes'};
    {'Professional2'}, {'house'};
    {'Professional2'}, {'vacuum_cleaner'};
    {'Professional3'}, {'bumps'};
    {'Professional3'}, {'house'};
    {'Professional3'}, {'mixer'};
    {'Professional3'}, {'tubes'};
    {'Professional4'}, {'hairdryer'};
    {'Professional4'}, {'vacuum_cleaner'};
    {'Professional5'}, {'house'};
    {'Professional5'}, {'wobble_surface'};
    {'Professional6'}, {'hairdryer'};    
    {'student1'}, {'house'};
    {'student1'}, {'shampoo_bottle'};
    {'student2'}, {'house'};
    {'student3'}, {'wobble_surface'};
    {'student4'}, {'house'};
    {'student5'}, {'house'};
    {'student5'}, {'wobble_surface'};
    {'student6'}, {'house'};
    {'student7'}, {'vacuum_cleaner'};
    {'student8'}, {'hairdryer'};
    {'student9'}, {'bumps'};
    {'student9'}, {'tubes'};
];


path_jerry = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\annotations\clean_all_jerry_both_passes';
path_dave_14 = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\annotations\dave-cleaned-pass1\cleaned-pass1';
path_dave_rest = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\annotations\dave-cleaned-others\cleaned-others';
path_yuan_rest = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\annotations\clean_rest_yuan\clean_rest_yuan';
path_original = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\full_objs';
folder_save = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v50rough_400\wrong_strokes_labelling\svgs_highlighted_errors';

if ~exist(folder_save)
    mkdir(folder_save);
end
% T_14 = table(Designer,Shape,Dave,Jerry);
% T_rest = table(Designer,Shape,Jerry,Yuan);

%% Read t 14 original files
num_strokes_file = zeros(size(list_designers_paper_14,1),1);
strokes_numbers_original = cell(size(list_designers_paper_14,1),1);
for i = 1:size(list_designers_paper_14,1)
    desinger{i} = list_designers_paper_14{i,1};
    object{i} = list_designers_paper_14{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough.obj'];
    strokes_numbers_original{i} = read_files(path_original, filename);
    num_strokes_file(i) = length(strokes_numbers_original{i});
end    

%% Read t 14 original files dave

num_strks_fl_dave = zeros(size(list_designers_paper_14,1),1);
strokes_numbers_dave = cell(size(list_designers_paper_14,1),1);
strokes_numbers_dave_removed = cell(size(list_designers_paper_14,1),1);
strokes_numbers_dave_kept = cell(size(list_designers_paper_14,1),1);
num_strks_dave_removed = zeros(size(list_designers_paper_14,1),1);
percentage_removed_dave = zeros(size(list_designers_paper_14,1),1);

for i = 1:size(list_designers_paper_14,1)
    desinger{i} = list_designers_paper_14{i,1};
    object{i} = list_designers_paper_14{i,2};
    
    filename = [desinger{i} '_' object{i} '_clean.obj'];
    strokes_numbers_dave{i} = read_files(path_dave_14, filename);
    num_strks_fl_dave(i) = length(strokes_numbers_dave{i});
    strokes_numbers_dave_removed{i} = setdiff(strokes_numbers_original{i}, strokes_numbers_dave{i});
    strokes_numbers_dave{i} = setdiff(strokes_numbers_original{i}, strokes_numbers_dave_removed{i});
    strokes_numbers_dave_kept{i} = setdiff(strokes_numbers_original{i}, strokes_numbers_dave_removed{i});
    
    num_strks_dave_removed(i) = length(strokes_numbers_dave_removed{i});
    percentage_removed_dave(i) = num_strks_dave_removed(i)/num_strokes_file(i)*100;   
    
end    
mean(percentage_removed_dave)
median(percentage_removed_dave)
%% Read t 14 original files jerry
num_strks_fl_jerry = zeros(size(list_designers_paper_14,1),1);
strokes_numbers_jerry = cell(size(list_designers_paper_14,1),1);
strokes_numbers_jerry_removed = cell(size(list_designers_paper_14,1),1);
percentage_removed_jerry = zeros(size(list_designers_paper_14,1),1);
num_strks_jerry_removed = zeros(size(list_designers_paper_14,1),1);
percentage_removed_14 = zeros(size(list_designers_paper_14,1),1);

for i = 1:size(list_designers_paper_14,1)
    desinger{i} = list_designers_paper_14{i,1};
    object{i} = list_designers_paper_14{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough.obj'];
    strokes_numbers_jerry{i} = read_files(path_jerry, filename);
    num_strks_fl_jerry(i) = length(strokes_numbers_jerry{i});
    strokes_numbers_jerry_removed{i} = setdiff(strokes_numbers_original{i}, strokes_numbers_jerry{i});
    
    num_strks_jerry_removed(i) = length(strokes_numbers_jerry_removed{i});
    percentage_removed_jerry(i) = num_strks_jerry_removed(i)/num_strokes_file(i)*100;    
    
    percentage_removed_14(i) = 0.5*(percentage_removed_jerry(i)+percentage_removed_dave(i));
end    
mean(percentage_removed_jerry)
mean(percentage_removed_14)

t14 = table(desinger', object', num_strokes_file, num_strks_dave_removed, percentage_removed_dave, num_strks_jerry_removed, percentage_removed_jerry, percentage_removed_14);

%% Read t supp original
num_strokes_file_supp = zeros(size(list_designers_supplemental,1),1);
strokes_numbers_original_supp = cell(size(list_designers_supplemental,1),1);
for i = 1:size(list_designers_supplemental,1)
    desinger{i} = list_designers_supplemental{i,1};
    object{i} = list_designers_supplemental{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough.obj'];
    strokes_numbers_original_supp{i} = read_files(path_original, filename);
    num_strokes_file_supp(i) = length(strokes_numbers_original_supp{i});
end    

%% Read t supp jerry
num_strks_fl_jerry_supp = zeros(size(list_designers_supplemental,1),1);
strokes_numbers_jerry_supp = cell(size(list_designers_supplemental,1),1);
strokes_numbers_jerry_removed_supp = cell(size(list_designers_supplemental,1),1);

n_removed_jerry_supp = zeros(size(list_designers_supplemental,1),1);
percentage_removed_jerry_supp = zeros(size(list_designers_supplemental,1),1);

for i = 1:size(list_designers_supplemental,1)
    desinger{i} = list_designers_supplemental{i,1};
    object{i} = list_designers_supplemental{i,2};
    
    filename = [desinger{i} '_' object{i} '_rough.obj'];
    strokes_numbers_jerry_supp{i} = read_files(path_jerry, filename);
    num_strks_fl_jerry_supp(i) = length(strokes_numbers_jerry_supp{i});
    strokes_numbers_jerry_removed_supp{i} = setdiff(strokes_numbers_original_supp{i}, strokes_numbers_jerry_supp{i});
    n_removed_jerry_supp(i) = length(strokes_numbers_jerry_removed_supp{i});
    percentage_removed_jerry_supp(i) = n_removed_jerry_supp(i)/num_strokes_file_supp(i)*100;    
        
end    
fprintf('Mean_supp_additional_jerry = %.2f\n', mean(percentage_removed_jerry_supp))

%% Read t supp yuan
num_strks_fl_yuan_supp = zeros(size(list_designers_supplemental,1),1);
strokes_numbers_yuan_supp = cell(size(list_designers_supplemental,1),1);
strokes_numbers_yuan_removed_supp = cell(size(list_designers_supplemental,1),1);
percentage_removed_yuan_supp = zeros(size(list_designers_supplemental,1),1);
n_removed_yuan_supp =  zeros(size(list_designers_supplemental,1),1);
percentage_removed_supp = zeros(size(list_designers_supplemental,1),1);

for i = 1:size(list_designers_supplemental,1)
    desinger{i} = list_designers_supplemental{i,1};
    object{i} = list_designers_supplemental{i,2};
    
    filename = [desinger{i} '_' object{i} '_clean.obj'];
    strokes_numbers_yuan_supp{i} = read_files(path_yuan_rest, filename);
    num_strks_fl_yuan_supp(i) = length(strokes_numbers_yuan_supp{i});
    strokes_numbers_yuan_removed_supp{i} = setdiff(strokes_numbers_original_supp{i}, strokes_numbers_yuan_supp{i});
    n_removed_yuan_supp(i) = length(strokes_numbers_yuan_removed_supp{i});
    
    percentage_removed_yuan_supp(i) = length(strokes_numbers_yuan_removed_supp{i})/num_strokes_file_supp(i)*100;
    
    
    percentage_removed_supp(i) = 0.5*(percentage_removed_yuan_supp(i)+percentage_removed_jerry_supp(i));        
end    
fprintf('Mean_supp_additional_yuan = %.2f\n', mean(percentage_removed_yuan_supp));
fprintf('Mean_supp_additional = %.2f\n', mean(percentage_removed_supp));

percentage_removed = [percentage_removed_14; percentage_removed_supp];

fprintf('Mean_supp = %.2f\n', mean([percentage_removed_14; percentage_removed_supp]));
fprintf('Median_supp = %.2f\n', median([percentage_removed_14; percentage_removed_supp]));

t_supp = table(desinger', object', num_strokes_file_supp, n_removed_jerry_supp, percentage_removed_jerry_supp, n_removed_yuan_supp, percentage_removed_yuan_supp, percentage_removed_supp);

%% 
% analysis_assigned_last;

%%



