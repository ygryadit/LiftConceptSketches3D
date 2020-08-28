function doAgressiveAssignement = checkNumConfigurationsAnyStroke(strokes_topology,  candidate_lines, inds_non_assgnd_dpndnt_strks)

doAgressiveAssignement = false;

global thr_max_num_config;
try
    vals = cat(2, candidate_lines(:).configurations);
catch e
    rethrow(e);
end
if length(vals) > thr_max_num_config
    doAgressiveAssignement = true;
    return;
end

for li = inds_non_assgnd_dpndnt_strks
    if isfield(strokes_topology(li), 'candidate_lines') & ~isempty(strokes_topology(li).candidate_lines)
        num_config = 0;
        for j = 1:length(strokes_topology(li).candidate_lines)
            num_config = num_config + length(strokes_topology(li).candidate_lines(j).configurations);
        end
  
        if length(vals) > thr_max_num_config
            doAgressiveAssignement = true;
            return;
        end
    end
end




end