function p_plane = checkPlanarityGivenConfiguration(candidate_line_dir, ... the stoke which is being updated
                                            configuration, ... the configuratiosn in which the added intersections is
                                            intersections,...
                                            cur_stroke_ind,... the stoke which is being updated
                                            strokes_topology,...
                                            ind_inter)

%% Planarity with lines that have been assigned
p_plane = 0;
for i = 1:length(configuration.inds_intrsctns__assigned)
    ind_intrsctn = configuration.inds_intrsctns__assigned(i);
    
    ind_intrsctng_strk = setdiff(intersections(ind_intrsctn).strokes_indices, ...
                                  cur_stroke_ind);

        
%     try
        
    if isfield(strokes_topology, 'planes_normals') & ~isempty(strokes_topology(ind_intrsctng_strk).planes_normals)
        normals = selectNormalsExcludeCurStroke(...
            strokes_topology(ind_intrsctng_strk).planes_normals, ...
            cur_stroke_ind);
        
        p_plane_ = computeProbabilityLieInPlane(...
                        normals,...
                        candidate_line_dir);
    else
        p_plane_ = 0;
    end
%     catch e
%         rethrow(e)
%     end
    p_plane = max(p_plane, p_plane_);
end




%% Planarity with lines that have not been assigned
if ~exist('ind_inter', 'var')
    ind_inter = [];
end



for i = 1:length(configuration.inds_intrsctns__mult_cnddts)
   ind_intrsctn = configuration.inds_intrsctns__mult_cnddts(i);
   if ind_intrsctn  == ind_inter
       continue;
   end
   ind_intrsctng_strk = configuration.inds_jnts_strks(i);

       
   ind_intrsctn_vrsn = configuration.inds_intrsctns__mult_cnddts_ind(i);
   
   % Find which configurations that correspond
   % Go over all the configurations and evaluate planaarity
   try
       mhs = find(intersections(ind_intrsctn).strokes_indices == ind_intrsctng_strk);
       inds_cnd_lns = intersections(ind_intrsctn).cnddts3D(ind_intrsctn_vrsn).cnddt_lns{mhs};
   catch e
       rethrow(e);
   end
   
   for ii = 1:length(inds_cnd_lns)
      ind_cnddt_ln = inds_cnd_lns(ii);
      candidate_line = strokes_topology(ind_intrsctng_strk).candidate_lines(ind_cnddt_ln);
      try
        inds_cnfgrtns = intersections(ind_intrsctn).cnddts3D(ind_intrsctn_vrsn).cnfgrtns{mhs}{ii};
      catch e
          rethrow e
      end
      for iii = 1:length(inds_cnfgrtns)
          try
              if ~isempty(candidate_line.configurations(inds_cnfgrtns(iii)).planes_normals)
              normals = selectNormalsExcludeCurStroke(...
                            candidate_line.configurations(inds_cnfgrtns(iii)).planes_normals, ...
                            cur_stroke_ind);
        
              p_plane_ = computeProbabilityLieInPlane(...
                                normals,...
                                candidate_line_dir);
              else
                  p_plane_ = 0;
              end
          catch e
              rethrow(e)
          end
      end
      p_plane = max(p_plane, p_plane_);
   end
   
   
   
end




end

