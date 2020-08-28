% Description:
% Once all the intersections and comfigurations are filled in, upate the
% directional vetors and normals for newly create dintersections.
% 
% It is doen after, since at the moment of creation not all configurations
% for the curren stroke are filled in.

function strokes_topology = updateDirectionalPriorsNormals( ...
                                                        strokes_topology,...
                                                        intersections,...
                                                        candidate_lines,...
                                                        ind_strk_nw)
                                                    
% Go over canidate lines and find the strokes and their respective
% configurations where the current stroke appears;
% For each such configuration recompute the normals and update directional
% prior if the line stroke is not the stroke towards vanishing point.

for i = 1:length(candidate_lines)
    candidate_line_strk_nw = candidate_lines(i);
    
    for j =1:length(candidate_line_strk_nw.configurations)
        configuration_strk_nw = candidate_line_strk_nw.configurations(j);
        
        inds_intrsctns__mult_cnddts = ...
            cat(1,configuration_strk_nw.inds_intrsctns__mult_cnddts);

        inds_intrsctns__mult_cnddts_ind = ...
            cat(1,configuration_strk_nw.inds_intrsctns__mult_cnddts_ind);

        inds_jnts_strks = ...
            cat(1,configuration_strk_nw.inds_jnts_strks);

        for ii = 1:length(inds_intrsctns__mult_cnddts)
            % Find line and configurations to update:

            ind_intrsctn = inds_intrsctns__mult_cnddts(ii);
            ind_intrsctn_vrsn = inds_intrsctns__mult_cnddts_ind(ii);
            ind_intrctng_strk = inds_jnts_strks(ii);

            mhs = find( ...
                intersections(ind_intrsctn).strokes_indices == ...
                ind_intrctng_strk);

            intrctns_cnddts3D = intersections(ind_intrsctn).cnddts3D(ind_intrsctn_vrsn);
            cnd_lns = intrctns_cnddts3D.cnddt_lns{mhs};
            
            for iii = 1:length(cnd_lns)
               cnd_ln = cnd_lns(iii);
               
               if strokes_topology(ind_intrctng_strk).line_group == 4
                   %Update line's configurations:
                    strokes_topology(ind_intrctng_strk).candidate_lines(cnd_ln) = ...
                        updateDirectionalCost(strokes_topology(ind_intrctng_strk).candidate_lines(cnd_ln), ...
                                              candidate_line_strk_nw,...
                                              configuration_strk_nw,...
                                              intrctns_cnddts3D,...
                                              ind_intrctng_strk,...
                                              mhs,...
                                              iii);
               end
               
               %Update normals:
               try
               strokes_topology(ind_intrctng_strk).candidate_lines(cnd_ln) = ...
                   updateNormals(strokes_topology(ind_intrctng_strk).candidate_lines(cnd_ln), ...
                          candidate_line_strk_nw,...
                          intrctns_cnddts3D,...
                          mhs,...
                          iii,...
                          ind_strk_nw);
               catch e
                   fprintf('ind_intrsctn = %d\n', ind_intrsctn);
                   rethrow(e);
               end
                
            
            end
        
        end
        
    end

end


end


function candidate_line_strk_updt = updateNormals(candidate_line_strk_updt, ...
                          candidate_line_strk_nw,...
                          intrctns_cnddts3D,...
                          mhs,...
                          iii,...
                          ind_strk_nw)

 normal = cross(candidate_line_strk_updt.dir,...
                     candidate_line_strk_nw.dir);
               
if norm(normal) < 0.5
   return;
end

% Add normal to each of the configurations:
try
inds_cnfgrtns = intrctns_cnddts3D.cnfgrtns{mhs}{iii};
catch e
    rethrow(e);
    
end
    
for ki = 1:length(inds_cnfgrtns )
    configuration = candidate_line_strk_updt.configurations(inds_cnfgrtns(ki));
    configuration.planes_normals(end+1,:) = [normal ind_strk_nw];
    candidate_line_strk_updt.configurations(inds_cnfgrtns(ki)) = configuration;
end

end

function candidate_line_strk_updt = ...
    updateDirectionalCost(candidate_line_strk_updt, ...
                          candidate_line_strk_nw,...
                          configuration_strk_nw,...
                          intrctns_cnddts3D,...
                          ind_stroke_update,...
                          mhs,...
                          iii)

p_ortho = computeScoreOrthogonality(candidate_line_strk_updt.dir, ...
                                    candidate_line_strk_nw.dir);
p_tangent = computeScoreTangential(candidate_line_strk_updt.dir,...
                                    candidate_line_strk_nw.dir);  

inds_cnfgrtns = intrctns_cnddts3D.cnfgrtns{mhs}{iii};

for ki = 1:length(inds_cnfgrtns )
   configuration = candidate_line_strk_updt.configurations(inds_cnfgrtns(ki));
   
  if  ~isempty(configuration_strk_nw.planes_normals)
      normals = selectNormalsExcludeCurStroke(...
                configuration_strk_nw.planes_normals, ...
                ind_stroke_update);

      p_plane = computeProbabilityLieInPlane(...
                            normals,...
                            candidate_line_strk_updt.dir);
  else
      p_plane = 0;
  end

  p_coverage   = configuration.p_coverage;
  p_directional   = configuration.p_directional;
  candidate_line_strk_updt.configurations(inds_cnfgrtns(ki)).p_directional = ...
      max([p_directional, p_ortho, p_tangent, p_plane]);

  p_full = costComposite( p_coverage, ...
             p_directional,...
             4);
  candidate_line_strk_updt.configurations(inds_cnfgrtns(ki)).p_full = ...
      p_full;
end

end