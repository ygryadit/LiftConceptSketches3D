function [inds_non_assgnd_dpndnt_strks, visited_strks] = ...
    findNonAssgndDpndntStrksSimplified(cur_stroke,...                             
                             strokes_topology,...
                             intersections,...
                             visited_strks)

if ~isfield(cur_stroke, 'inds_intrsctns_eval_mltpl_cnddts')                         
    visited_strks = [visited_strks cur_stroke.ind];
    inds_non_assgnd_dpndnt_strks = [];
    return;
end
strks_pairs = [intersections(cur_stroke.inds_intrsctns_eval_mltpl_cnddts).strokes_indices];
inds_non_assgnd_dpndnt_strks = setdiff(strks_pairs(:), cur_stroke.ind);

inds_non_assgnd_dpndnt_strks = unique(inds_non_assgnd_dpndnt_strks);
mask = ~[strokes_topology(inds_non_assgnd_dpndnt_strks).depth_assigned];
inds_non_assgnd_dpndnt_strks = inds_non_assgnd_dpndnt_strks(mask);



non_visited_strks = setdiff(inds_non_assgnd_dpndnt_strks, visited_strks);

while ~isempty(non_visited_strks)
    
    dependent_stroke = strokes_topology(non_visited_strks(1));
    dependent_stroke.ind = non_visited_strks(1);
%     fprintf('Checking stroke %d\n', dependent_stroke.ind);
    visited_strks = [visited_strks dependent_stroke.ind];
    
    UP_TO_LAST = true;
    [dependent_stroke.inds_intrsctns_eval,...
     dependent_stroke.inds_intrsctns_eval_actv,...
     dependent_stroke.inds_intrsctns_eval_mltpl_cnddts,...
     dependent_stroke.inds_intrsctng_strks_eval,...
     dependent_stroke.inds_intrsctng_strks_eval_actv,...
     dependent_stroke.inds_intrsctng_strks_eval_mltpl_cnddts] = ...
        returnIndicesNodesTypes(dependent_stroke, ...
                            cat(1, strokes_topology(:).depth_assigned),...
                                        intersections,...
                                        UP_TO_LAST);
                                    
    [inds_non_assgnd_dpndnt_strks_, visited_strks_] = ...
        findNonAssgndDpndntStrksSimplified(dependent_stroke,...                             
                                 strokes_topology,...
                                 intersections,...
                                 visited_strks);
    
    inds_non_assgnd_dpndnt_strks = unique([inds_non_assgnd_dpndnt_strks; inds_non_assgnd_dpndnt_strks_]);
    visited_strks = unique([visited_strks visited_strks_]);
    
    non_visited_strks = setdiff(inds_non_assgnd_dpndnt_strks, visited_strks);
end

end
% 
% function inds_non_assgnd_dpndnt_strks = ...
%         findNonAssgndDpndntStrksReq(strokes_topology,...
%                                     intersections,...
%                                     confgrtns, ...
%                                     strokes_visited)
%     inds_non_assgnd_dpndnt_strks = [];
%     for iclc = 1:length(confgrtns)
%         confgrtn = confgrtns(iclc);
%         [inds_jnts_strks, ia] = setdiff(confgrtn.inds_jnts_strks, ...
%                                                strokes_visited);
%         
%         inds_intrsctns__mult_cnddts = confgrtn.inds_intrsctns__mult_cnddts(ia);
%         inds_intrsctns__mult_cnddts_ind = confgrtn.inds_intrsctns__mult_cnddts_ind(ia);
%         
%         
%         for j = 1:length(inds_jnts_strks)
%            strk = inds_jnts_strks(j);
%            
%            
%            intrsctn = intersections(inds_intrsctns__mult_cnddts(j));
%            ind_intrsctng_strk = find(intrsctn.strokes_indices ~= strokes_visited(end));
%            
%            inds_lnd = intrsctn.cnddts3D(inds_intrsctns__mult_cnddts_ind(j)).cnddt_lns{ind_intrsctng_strk};
%            
%            for i = 1:length(inds_lnd )
%            
%                candidate_line = strokes_topology(strk).candidate_lines(...
%                                         inds_lnd(i));
% 
%                configurations =  candidate_line.configurations(...
%                    intrsctn.cnddts3D(inds_intrsctns__mult_cnddts_ind(j)).cnfgrtns{ind_intrsctng_strk}{i});
% 
% 
%                vals = findNonAssgndDpndntStrksReq(strokes_topology,...
%                                     intersections,...
%                                     configurations, ...
%                                     [strokes_visited strk]);  
%                                 
%                inds_non_assgnd_dpndnt_strks = [inds_non_assgnd_dpndnt_strks,...
%                    vals];
%            end
%         end
%     end
% end