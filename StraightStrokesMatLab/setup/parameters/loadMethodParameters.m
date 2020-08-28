
%% Parameters
% global confidence_threshold;
% confidence_threshold = 0.2; %0.25

global confidence_threshold;
confidence_threshold = 0.2;

global confidence_threshold2;
confidence_threshold2 = 0.1;

global confidence_step;
confidence_step = 0.75;

global confidence_threshold_min;
confidence_threshold_min = 0.01;


global max_cost_threshold;
% max_cost_threshold = 0.65; %(0.4 + 0.9)/2
% max_cost_threshold = 0.75; 
% max_cost_threshold = 0.6; 
max_cost_threshold = 0.75; 

global max_cost_absolute;
max_cost_absolute = 0.98; %used to be 95
% max_cost_threshold = 0.85; 

global stokes_limit;
stokes_limit = Inf;

global dist_endpoint;
% dist_endpoint = 0.15;
dist_endpoint = 0.25;

global thr_common; % threshold to merge the lines (the precentage of common length from the length of the shortest)
thr_common = 0.85;

global sigma_costs;
sigma_costs = cos(pi/2 - pi/36);

global sigma_aligned;
sigma_aligned = sind(0.75);
%% Speedup

global thr_max_num_lines;
% thr_max_num_lines = 200;
thr_max_num_lines = 200;

global thr_max_num_config;
thr_max_num_config = 500;

%% Variables
global last_added_stroke;
last_added_stroke = 1;

global CameraView;
CameraView = [-64.1805,   25.6400];

global frame_count;
frame_count = 0;