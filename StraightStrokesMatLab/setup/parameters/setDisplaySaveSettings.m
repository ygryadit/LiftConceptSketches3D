global DISPLAY_INFO
DISPLAY_INFO = false;

global DISPLAY_STATUS
DISPLAY_STATUS = false;

global SHOW_FIGS
SHOW_FIGS = false;


global SHOW_FIGS_PREPROCESS
SHOW_FIGS_PREPROCESS = false;

global DO_BUILD_GRAPH;
DO_BUILD_GRAPH = false;

global DEBUG;
DEBUG = false;

if ~DEBUG
   warning off;
end



%% Save
global SAVE_PROCESS;
SAVE_PROCESS = false;

global SAVE_SVGS;
SAVE_SVGS = true;