% Example usage:
% launchOpenSketch('D:\Projects\SketchItProject\CollectedSketches\SKETCHES_PUBLISH\sketches_json_first_viewpoint\',...
% 'student8',...
% 'house',...
% 'Z:\WiresProject\Code\Matlab\reconstruct_wires\output_new')

function launchOpenSketch_cluster_a(folder_designer_,...
                        designer_,...
                        object_name_,...
                        folder_out,...
                        view_, ... %'view1' or 'view2'
                        datatset_name_) %'OpenSketch'

if ~exist('view_', 'var')
    view_ = 'view1';
end

if ~exist('datatset_name_', 'var')
    datatset_name_ = 'OpenSketch';
end
%% Globals
global datatset_name;
datatset_name = datatset_name_;
    
if strcmp(datatset_name, 'OpenSketch')
    global folder_designer;
    folder_designer = folder_designer_;
    global designer;
    designer = designer_;
    global object_name;
    object_name = object_name_;
    global view;
    view = view_;               
end

global folder_save;        
folder_save = folder_out;
global folder_save_imgs;
global filepath_sketch_json;
global filepath_sketch_img;
global folderVP;

%% Method parameters:
loadMethodParameters;


global thr_max_num_lines;
thr_max_num_lines = 400;


%% Visualization parameters:
setDisplaySaveSettings;

%% Behaviour:
setMethodOptions;

%% Logfile
global fid;
% fid = fopen(fullfile(folder_save, 'LogFile.txt'), 'a');
% if fid == -1
%   error('Cannot open log file.');
% end
% fprintf(fid, 'Ready \n');
fid = 1;
% fid 1 to output to matlab screen.
fprintf(fid, 'Designer %s, object %s \n', designer, object_name);


%% Folders:
folder_save = fullfile(folder_out, 'gen_cand_planes', sprintf('_%d', thr_max_num_lines));

[folder_save, ...
 folder_save_imgs, ...
 filepath_sketch_json,...
 filepath_sketch_img] = setupFolderPaths(folder_save);


folderVP = fullfile(folder_save, 'EstimatedVP');
if ~exist(folderVP, 'dir')
    mkdir(folderVP);
end 

%% Run
timerVal1 = tic;
global GENERATE_CANDIDATE_PLANES;
GENERATE_CANDIDATE_PLANES = true;



try
    roughSketch3DInference();
    fprintf(fid, 'Success: designer %s, object %s \n', designer, object_name);
    ellapsed_time = toc(timerVal1);
    fprintf(fid,'\t\t Sketch time %.3f\n', ellapsed_time);
    save(fullfile(folder_save, 'preformance.mat'), 'ellapsed_time');

catch e
    getReport(e, 'extended')
end                  

% fclose all;


end