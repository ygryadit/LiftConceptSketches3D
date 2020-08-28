function [inds_non_assgnd_dpndnt_strks] = ...
    findNonAssgndDpndntStrks(cur_stroke,...
                             candidate_lines,...
                             strokes_topology,...
                             intersections)

inds_non_assgnd_dpndnt_strks = [];
strokes_visited = cur_stroke.ind;
for icl = 1:length(candidate_lines)
    
    confgrtns = candidate_lines(icl).configurations;
    
    for iclc = 1:length(confgrtns)
        confgrtn = confgrtns(iclc);
        
        for iimc = 1:length(confgrtn.inds_intrsctns__mult_cnddts)
            
            inds_non_assgnd_dpndnt_strks = [inds_non_assgnd_dpndnt_strks, ...
                                            confgrtn.inds_jnts_strks];
            
            for j = 1:length(confgrtn.inds_jnts_strks)
                strk = confgrtn.inds_jnts_strks(j);
                candidate_line = strokes_topology(strk).candidate_lines(...
                                    confgrtn.inds_intrsctns__mult_cnddts_ind(j));
                
%                 vals = findNonAssgndDpndntStrksReq(strokes_topology,...
%                                     intersections,...
%                                     candidate_line.configurations, ...
%                                     [strokes_visited strk]);
                inds_non_assgnd_dpndnt_strks = [inds_non_assgnd_dpndnt_strks, ...
                    candidate_line.list_jnt_strks];
            end
            
            
        end
        
    end
end

inds_non_assgnd_dpndnt_strks = unique(inds_non_assgnd_dpndnt_strks);
mask = ~[strokes_topology(inds_non_assgnd_dpndnt_strks).depth_assigned];
inds_non_assgnd_dpndnt_strks = inds_non_assgnd_dpndnt_strks(mask);

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