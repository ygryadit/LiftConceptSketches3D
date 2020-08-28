global DISPLAY_INFO
DISPLAY_INFO = true;

global DISPLAY_STATUS
DISPLAY_STATUS = false;

global SHOW_FIGS
SHOW_FIGS = true;


global SHOW_FIGS_PREPROCESS
SHOW_FIGS_PREPROCESS = true;

global DO_BUILD_GRAPH;
DO_BUILD_GRAPH = false;

global DEBUG;
DEBUG = true;

if ~DEBUG
   warning off;
end



%% Save
global SAVE_PROCESS;
SAVE_PROCESS = false;

global SAVE_SVGS;
SAVE_SVGS = true;