function doMerge = checkMerge(begin, end_coordinate,strokes_topology,ind_strk_1, ind_strk_2, dist, areContinious_)%
% l_common, length_s1, length_s2,strokes_topology,ind_strk_1, ind_strk_2, dist)

    % Length of the common region:
   global thr_common; % threshold to merge the lines (the precentage od common length from the length of the shortest)
   thr_common = 0.75; 
   l_common = norm(end_coordinate - begin);
   
   
%    doMerge = (dist < 0.1*strokes_topology(ind_strk_1).length2DFull) & ...
%              (dist < 0.1*strokes_topology(ind_strk_2).length2DFull);
    
%     doMerge = (dist < strokes_topology(ind_strk_1).accuracy_radius) & ...
%               (dist < strokes_topology(ind_strk_2).accuracy_radius);
    
%     doMerge = (dist < (strokes_topology(ind_strk_2).accuracy_radius + strokes_topology(ind_strk_1).accuracy_radius));

%     cond_legth_1 = l_common/ length_s1  >= thr_common ;
%     cond_legth_2 = l_common/ length_s2  >= thr_common ;
%     
   cond_legth_1 = l_common/ strokes_topology(ind_strk_1).length2DPrimitive  >= thr_common ;
    cond_legth_2 = l_common/ strokes_topology(ind_strk_2).length2DPrimitive  >= thr_common ;
%     r = strokes_topology(ind_strk_2).accuracy_radius;
    r = max([strokes_topology(ind_strk_2).accuracy_radius strokes_topology(ind_strk_1).accuracy_radius]);
    
    areClose = (dist < r);
    
    
    areSimilarLengths = (cond_legth_1 & cond_legth_2);
    areStronglyOverlaping = (cond_legth_1 | cond_legth_2);
        
    areContinious = (abs(ind_strk_2 - ind_strk_1) <= 2) | areContinious_;
    
    %Can add the check that there is no earlier intersecting stoke
    
    
    
    areCloseInTime = (abs(ind_strk_2 - ind_strk_1) <= 5);
    
    doMerge = areClose & (....
                    areContinious | ...
                    areSimilarLengths | ...
                    (areCloseInTime & areStronglyOverlaping) ...
                );




% 
%     doMerge = ((abs(ind_strk_2 - ind_strk_1) <= 5) | areSimilarLengths);
%     
%     doMerge =  areClose & doMerge;
%           
%     doMerge = doMerge & areStronglyOverlaping ;
end