% When the line can not be assigned we same the candidate lines, and add
% the paired information into candidate intersections.

function intersections = assignCandidateIntersectionSecondStroke(candidate_lines, intersections, cur_stroke_ind)
    
for cli = 1:length(candidate_lines)
   for ccj = 1:length(candidate_lines(cli).configurations)
       configuration = candidate_lines(cli).configurations(ccj);
       
       for iim = 1:length(configuration.inds_intrsctns__mult_cnddts)
           ind_intersection  = configuration.inds_intrsctns__mult_cnddts(iim);
           ind_candidate  = configuration.inds_intrsctns__mult_cnddts_ind(iim);
           
           mhs = find(cur_stroke_ind == intersections(ind_intersection).strokes_indices);
           
           
           ind_member = find(ismember(intersections(ind_intersection).cnddts3D(ind_candidate).cnddt_lns{mhs}, cli));
      
           if ~isempty(ind_member)
               %add configration
               intersections(ind_intersection).cnddts3D(ind_candidate).cnfgrtns{mhs}{ind_member}(end+1) = ccj;
               intersections(ind_intersection).cnddts3D(ind_candidate).cnfgrtns{mhs}{ind_member} = ...
                   unique(intersections(ind_intersection).cnddts3D(ind_candidate).cnfgrtns{mhs}{ind_member});
           else
               %create new
           
               ind_cli = length(intersections(ind_intersection).cnddts3D(ind_candidate).cnddt_lns{mhs})+1;
               intersections(ind_intersection).cnddts3D(ind_candidate).cnddt_lns{mhs}(ind_cli) = cli;

               if isempty(intersections(ind_intersection).cnddts3D(ind_candidate).cnfgrtns{mhs})
                    intersections(ind_intersection).cnddts3D(ind_candidate).cnfgrtns{mhs}{1}(1) = ccj;           
               else
                   if length(intersections(ind_intersection).cnddts3D(ind_candidate).cnfgrtns{mhs}) < ind_cli || ...
                         isempty(intersections(ind_intersection).cnddts3D(ind_candidate).cnfgrtns{mhs}{ind_cli})
                       intersections(ind_intersection).cnddts3D(ind_candidate).cnfgrtns{mhs}{ind_cli} = [];
                   end
                   intersections(ind_intersection).cnddts3D(ind_candidate).cnfgrtns{mhs}{ind_cli}(end+1) = ccj;
               end
           end
           %fprintf('ind_intersection %d, ind_candidate %d, mhs %d, ind_cli %d\n', ind_intersection, ind_candidate, mhs, ind_cli)
       end
   end
end

end