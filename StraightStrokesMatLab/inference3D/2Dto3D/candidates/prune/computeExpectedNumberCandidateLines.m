function val = computeExpectedNumberCandidateLines(cur_stroke, strokes_topology, intersections)

k = length(cur_stroke.inds_intrsctns_eval_actv);

val = k*(k-1)/2;

strks_pairs = [intersections(cur_stroke.inds_intrsctns_eval_mltpl_cnddts).strokes_indices];
Inds_strks_mltpl_cnddts = setdiff(strks_pairs(:), cur_stroke.ind);

num_strks = length(Inds_strks_mltpl_cnddts);
M = zeros(1,num_strks);
for i = 1:num_strks
    ind_strk = Inds_strks_mltpl_cnddts(i);
    if isfield(strokes_topology(ind_strk), 'candidate_lines')
        M(i) = length(strokes_topology(ind_strk).candidate_lines);
    else
        M(i) = 0;
    end
end
 
val = val + sum(k*M);


for i = 1:num_strks
   inds = [1:(i-1) (i+1):num_strks];
   val = val + M(i)*sum(M(inds));
end 
end