function candidate_line = addIntersectionsConfigCandidateLine(...
                                candidate_line,...
                                inds_intrsctns,...
                                inds_intrsctns__assigned,...
                                inds_intrsctns__mult_cnddts,...
                                inds_intrsctns__mult_cnddts_ind,...
                                inds_jnts_strks,...                                
                                p_intrsctns_dists,...
                                p_directional,...
                                p_coverage,...
                                p_full,...
                                planes_normals)
                            
    if ~isfield(candidate_line, 'configurations')
        ind_cnfgrtn = 1;
    else
        ind_cnfgrtn = length(candidate_line.configurations) +1;
    end
    
    candidate_line.configurations(ind_cnfgrtn).inds_intrsctns = inds_intrsctns;
    candidate_line.configurations(ind_cnfgrtn).inds_intrsctns__assigned = inds_intrsctns__assigned;
    candidate_line.configurations(ind_cnfgrtn).inds_intrsctns__mult_cnddts = inds_intrsctns__mult_cnddts;
    candidate_line.configurations(ind_cnfgrtn).inds_intrsctns__mult_cnddts_ind = inds_intrsctns__mult_cnddts_ind; %index in the intersections stuctre
%     candidate_line.configurations(ind_cnfgrtn).inds_jnts_strks_cnddte_lns = inds_jnts_strks_cnddte_lns;
 
    candidate_line.configurations(ind_cnfgrtn).inds_jnts_strks = inds_jnts_strks;
    
    candidate_line.configurations(ind_cnfgrtn).p_intrsctns_dists = p_intrsctns_dists;
    candidate_line.configurations(ind_cnfgrtn).p_directional = p_directional;
    candidate_line.configurations(ind_cnfgrtn).p_coverage = p_coverage;
    candidate_line.configurations(ind_cnfgrtn).p_full = p_full;
    if ~exist('planes_normals', 'var')
        planes_normals =[];
    end
    
    
    candidate_line.configurations(ind_cnfgrtn).p_full_joint = 0;
    candidate_line.configurations(ind_cnfgrtn).pfj_assgnd =[];
    candidate_line.configurations(ind_cnfgrtn).num_joint_strokes =[];     
    candidate_line.configurations(ind_cnfgrtn).list_jnt_strks =[];
    candidate_line.configurations(ind_cnfgrtn).planes_normals = planes_normals;
%     candidate_line.configurations(ind_cnfgrtn).p_joint = p_joint;
end