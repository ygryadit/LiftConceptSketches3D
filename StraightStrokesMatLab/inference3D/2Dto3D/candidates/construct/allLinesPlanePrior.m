% OUTPUT:
%     candidate_lines.origin
%     candidate_lines.dirs
%     candidate_lines.coordinates3D_prior //endpoints coordinates of
%           the lines proxy in case of the line.
% 
%     candidate_lines.p_directional //3D  
% 
%     candidate_lines.inds_jnts_strks_cnddte_lns
%     candidate_lines.inds_jnts_strks
%     candidate_lines.inds_intrsctns_jnts_strks
% 
% 
%     candidate_lines.configurations  
%           candidate_lines.configurations{i}.intersections
%           candidate_lines.configurations{i}.intersections__assigned
%           candidate_lines.configurations{i}.intersections__mult_cnddts
% 
%           candidate_lines.configurations{i}.p_intrsctns_dists
% 
%           candidate_lines.configurations{i}.p_coverage //2D  
%           candidate_lines.configurations{i}.p_full
% -------------------------------------------------------------------------
% Description:
% 
%   Genrates all the lines directly passing through each of the
%   intersection and lies in the given plane with a given projection.
    
function candidate_lines...
            = allLinesPlanePrior( ...
                    intersections,...
                    strokes_topology,...
                    cam_param,...
                    cur_stroke)
    
    %% First all the lines passing through the strokes with a selected 3D value:            
    inds_intrsctns_assigned = cur_stroke.inds_intrsctns_eval_actv;

    orth_vec = getDirectionVec(cur_stroke.ind_orth_ax);
    
    for i = 1:length(inds_intrsctns_assigned)
        
        indcs_intrsctns = inds_intrsctns_assigned(i);
        indcs_intrsctns__assigned = inds_intrsctns_assigned(i);
        indcs_intrsctns__mult_cnddts = [];%zeros(1,num_intrsctns_mult_hypoth)
        indcs_intrsctns__mult_cnddts_ind = [];
        inds_jnts_strks = [];                                
        
        
        p3D_in = intersections(inds_intrsctns_assigned(i)).coordinates3D;
        
        
        candidate_lines(i) = createOneCandidateLineGivenPlanePriorAndProjection(...
                                                    p3D_in,...
                                                    orth_vec,...
                                                    cur_stroke.primitive_geom,...
                                                    cam_param,...
                                                    indcs_intrsctns,...
                                                    indcs_intrsctns__assigned,...
                                                    indcs_intrsctns__mult_cnddts,...
                                                    indcs_intrsctns__mult_cnddts_ind,...
                                                    inds_jnts_strks);     
    end
    
    
    
    %% Lines originating from intersections with strokes with multiple possible solutions:
    
    global ESTIMATE_JOINTLY; 
    if ESTIMATE_JOINTLY       

        mask_prev_stokes_intersections  = true(size(cur_stroke.inds_intrsctns_eval_mltpl_cnddts));
        inds_intrsctns__mult_candidates = cur_stroke.inds_intrsctns_eval_mltpl_cnddts;
        
        inds_strokes__mult_candidates = cur_stroke.inds_intrsctng_strks_eval_mltpl_cnddts(mask_prev_stokes_intersections); 
        
        if exist('candidate_lines', 'var')
            lines_ind = length(candidate_lines);
        else
            lines_ind = 0;
        end
        
        for j = 1:length(inds_intrsctns__mult_candidates)
            
            ind_hypoth_inter = inds_intrsctns__mult_candidates(j);
            ind_hypoth_stroke = inds_strokes__mult_candidates(j);
            
            if ~isfield(strokes_topology(ind_hypoth_stroke), 'candidate_lines')
                continue;
            end

            % For each possible position of the 3D intersection create
            % candidate line:
            
            num_candidate_lines__ind_hypoth_stroke = length(strokes_topology(ind_hypoth_stroke).candidate_lines);
            for hj = 1:num_candidate_lines__ind_hypoth_stroke 
               
%                fprintf(2, 'It seems that the intersection in the line below shoulf be already computed');
         
               
               [inter_coord_3D, ~] = ...
                    opt3Dpos2DProj( cat(1,intersections(ind_hypoth_inter).coordinates2D), cam_param, ...
                                strokes_topology(ind_hypoth_stroke).candidate_lines(hj).coordinates3D_prior(1, 1:3), ...
                                strokes_topology(ind_hypoth_stroke).candidate_lines(hj).coordinates3D_prior(1, 4:6));
               
                lines_ind = lines_ind + 1;            
                
                indcs_intrsctns = ind_hypoth_inter;
                indcs_intrsctns__assigned = [];
                indcs_intrsctns__mult_cnddts = ind_hypoth_inter;
                indcs_intrsctns__mult_cnddts_ind = hj;
                inds_jnts_strks = ind_hypoth_stroke;                                
                
                candidate_lines(lines_ind) = ...
                    createOneCandidateLineGivenPlanePriorAndProjection(...
                                                    inter_coord_3D,...
                                                    orth_vec,...
                                                    cur_stroke.primitive_geom,...
                                                    cam_param,...
                                                    indcs_intrsctns,...
                                                    indcs_intrsctns__assigned,...
                                                    indcs_intrsctns__mult_cnddts,...
                                                    indcs_intrsctns__mult_cnddts_ind,...
                                                    inds_jnts_strks);
                                    
            end
        end
    end
    
    
    if ~exist('candidate_lines', 'var')
        candidate_lines = [];
    end
    
    
end
