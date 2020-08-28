global DIRECTIONAL_NON_AXIS;
DIRECTIONAL_NON_AXIS = true;

global USE_ONLY_LIKELY_INTERS;
USE_ONLY_LIKELY_INTERS = true;

global MERGE_LINES;
MERGE_LINES = true;

global ACCOUNT_INTER_CURVS;
ACCOUNT_INTER_CURVS = true;

global ESTIMATE_JOINTLY;
ESTIMATE_JOINTLY = true; 

% Do not assign depth if the stroke is not confident:
global DELAY_ASSIGN
DELAY_ASSIGN = true;

global ASSIGN_BEST;
ASSIGN_BEST = false;

%The parametet used in doNotAssignDepthValueToLine:
global ASSIGN_HIGH_SCORE;
ASSIGN_HIGH_SCORE = true;

global DO_PRUNE;
DO_PRUNE  = true;

global ORTHOGRAPHIC;
ORTHOGRAPHIC = false;

global USE_ADDITIONA_VPs;
USE_ADDITIONAL_VPs = false;

global USE_CURVES;
USE_CURVES = false;


global GENERATE_CANDIDATE_PLANES;
GENERATE_CANDIDATE_PLANES = false;

global USEMAXDIR;
USEMAXDIR = true;
 
global USE_SIMPLE_COST;
USE_SIMPLE_COST = false;
  
global USE_WEAK_REJECT;
USE_WEAK_REJECT = true;
    