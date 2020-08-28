function candidate_lines = initialiseCandidateLines(num_lines, ...
                                    inters_coord3D_begin,...
                                    inters_coord3D_end,...
                                    ind_inter_pairs,...
                                    inds_intrsctns__assigned,...
                                    inds_intrsctns__mult_cnddts,...
                                    inds_intrsctns__mult_cnddts_ind,...
                                    inds_jnts_strks,...
                                    cur_stroke,...
                                    cam_param)
                                
    dirs = inters_coord3D_begin - inters_coord3D_end;
    dirs = dirs./sqrt(sum(dirs.^2,2));                          
    
    % Set origins to the center points of the lines:
    origin = 0.5*(inters_coord3D_begin + inters_coord3D_end);
   
    for i = 1:num_lines
        
        candidate_lines(i) = initialiseOneCandidateLine(origin(i,:),...
                                                    dirs(i,:),...
                                                    cur_stroke.primitive_geom,...
                                                    cam_param.P);

        inds_intrsctns = ind_inter_pairs(i,:);
        
        
        
        if ~isempty(inds_intrsctns__assigned)
            inds_intrsctns__assigned_i = inds_intrsctns__assigned(i,:);
        else
            inds_intrsctns__assigned_i = [];
        end
            
        if ~isempty(inds_intrsctns__mult_cnddts)
            inds_intrsctns__mult_cnddts_i = inds_intrsctns__mult_cnddts(i,:);
        else
            inds_intrsctns__mult_cnddts_i = [];
        end
        
        if ~isempty(inds_intrsctns__mult_cnddts_ind)
            inds_intrsctns__mult_cnddts_ind_i = inds_intrsctns__mult_cnddts_ind(i,:);
        else
            inds_intrsctns__mult_cnddts_ind_i = [];
        end
        
        try
            if ~isempty(inds_jnts_strks)
                inds_jnts_strks_i = inds_jnts_strks{i};
            else
                inds_jnts_strks_i = [];
            end
        catch
            warning('test');
        end
%         inds_intrsctns__assigned = inds_intrsctns;
  
        p_intrsctns_dists = [1.0 1.0];
        p_directional = NaN;
        p_coverage = NaN;
        p_full = NaN;
        
        candidate_lines(i) = addIntersectionsConfigCandidateLine(...
                                candidate_lines(i),...
                                inds_intrsctns,...
                                inds_intrsctns__assigned_i,...
                                inds_intrsctns__mult_cnddts_i,...
                                inds_intrsctns__mult_cnddts_ind_i,...
                                inds_jnts_strks_i,...                                
                                p_intrsctns_dists,...
                                p_directional,...
                                p_coverage,...
                                p_full);         
    end
end