function [intrsctns_cmmn, ... % common intersection bewtween cnfigurtions.
          intrsctns_mlt_hpthss, ... % common intersection bewtween configurtions that are with lines with multiple candidates.
          strks_mlt_hpthss,... % strokes coorespondin to the intersections defined aboves
          inds_cnddt_intrsctn]... % common intersection version (ind) bewtween configurtions that are with lines with multiple candidates
            = findCommonSubsetOfIntersections(configurations, intersections, strokes_topology)
    
    intrsctns_cmmn = [];
    intrsctns_mlt_hpthss = [];
    strks_mlt_hpthss = [];
    inds_cnddt_intrsctn  = [];    
        
    %%     
    num_cnfgrtns = length(configurations);                 
    
    lines_strings = strings(num_cnfgrtns,1);

    for l_ind = 1:num_cnfgrtns
      if ~isempty(configurations(l_ind).inds_intrsctns__assigned)  
          lines_strings(l_ind) = ...
              sprintf('%d_', [configurations(l_ind).inds_intrsctns__assigned]);
      else
          lines_strings(l_ind) = '';
      end
      
      for j = 1:length(configurations(l_ind).inds_intrsctns__mult_cnddts)
          lines_strings(l_ind) =   sprintf('%s_%d_%d_', ...
                                        lines_strings(l_ind),...
                                        configurations(l_ind).inds_intrsctns__mult_cnddts(j), ...
                                        configurations(l_ind).inds_intrsctns__mult_cnddts_ind(j));
           
      end
      % Intersection indices and intersection versions.
    end

    [~, ind_unique, ~] = unique(lines_strings, ...
                                'rows', 'stable');
    
    % There is only one configuration:
    if length(ind_unique) == 1
         [intrsctns_cmmn, ... 
          intrsctns_mlt_hpthss, ... 
          strks_mlt_hpthss,... 
          inds_cnddt_intrsctn]... 
            = assignCommonSetsData(intersections, configurations, 1, strokes_topology);
       return;
    end
    
    % In case of multiple configurations find unique subset:
   
    for i = 1:num_cnfgrtns
        issubsetstr = true;
        for j = [1:i-1 i+1:num_cnfgrtns]
            issubsetstr  = issubsetstr & contains(lines_strings(j), lines_strings(i));
        end
        if issubsetstr
            break;
        end
    end
    
    if (i == num_cnfgrtns) && ~issubsetstr
        return;
    end
    
    [intrsctns_cmmn, ...
     intrsctns_mlt_hpthss, ... 
     strks_mlt_hpthss,... 
     inds_cnddt_intrsctn]... 
            = assignCommonSetsData(intersections, configurations, i, strokes_topology);
end


function  [intrsctns_cmmn, ... % common intersection bewtween cnfigurtions.
          intrsctns_mlt_hpthss, ... % common intersection bewtween configurtions that are with lines with multiple candidates.
          strks_mlt_hpthss,... % strokes coorespondin to the intersections defined aboves
          inds_cnddt_intrsctn]... % common intersection version (ind) bewtween configurtions that are with lines with multiple candidates
            = assignCommonSetsData(intersections, configurations, i, strokes_topology)

   intrsctns_cmmn = configurations(i).inds_intrsctns;

   mask_nodes_mult_hpthss = ismember(configurations(i).inds_intrsctns__mult_cnddts, intrsctns_cmmn);

   mask_not_activated = reshape(isnan(cat(1,intersections(configurations(i).inds_intrsctns__mult_cnddts).is_active)), 1,[]);
   
   mask_not_consdrd_intrsctns_multpl_cnddts = mask_nodes_mult_hpthss & mask_not_activated;
   
   inds_cnddt_intrsctn  = configurations(i).inds_intrsctns__mult_cnddts_ind(mask_not_consdrd_intrsctns_multpl_cnddts);
   strks_mlt_hpthss     = configurations(i).inds_jnts_strks(mask_not_consdrd_intrsctns_multpl_cnddts);    
   intrsctns_mlt_hpthss = configurations(i).inds_intrsctns__mult_cnddts(mask_not_consdrd_intrsctns_multpl_cnddts);
       
   
   inds_migrate = cat(1,strokes_topology(strks_mlt_hpthss).depth_assigned);
   
   intrsctns_cmmn = [intrsctns_cmmn intrsctns_mlt_hpthss(inds_migrate)];
   intrsctns_mlt_hpthss = intrsctns_mlt_hpthss(~inds_migrate);
   strks_mlt_hpthss = strks_mlt_hpthss(~inds_migrate);
   inds_cnddt_intrsctn = inds_cnddt_intrsctn(~inds_migrate);
end