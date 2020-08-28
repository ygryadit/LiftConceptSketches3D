function [candidate_lines, p_full_cnddt_lns, p_joint_cnddt_lns, p_joint, strokes_topology] = ...
                    computeJointCostsV2(strokes_topology,...
                                      intersections,...
                                      candidate_lines,...
                                      pDirectional,...
                                      cur_stroke)
    
    global ESTIMATE_JOINTLY;

    hi = [];

    debug_joint_cost = false;
    
    
    %% All other cost:
    num_lines = length(candidate_lines);
    
    p_joint_cnddt_lns = zeros(1,num_lines);
    p_joint = cell(num_lines,1);
    p_full_cnddt_lns = zeros(1,num_lines);

    for i = 1:num_lines         
        for j = 1:length(candidate_lines(i).configurations)
    
            if debug_joint_cost
                fig_num_3D = figure(10);
                plotStrokesTopolgyIntersectionsTypes(...
                                strokes_topology,...
                                cur_stroke,...
                                intersections,...
                                fig_num_3D)
                            
                plot3( candidate_lines(i).coordinates3D_prior(1,[1,4]),...
                         candidate_lines(i).coordinates3D_prior(1,[2,5]),...
                         candidate_lines(i).coordinates3D_prior(1,[3,6]),...
                         ':', 'LineWidth', 2);
                for h = 1:length(candidate_lines(i).configurations(j).inds_intrsctns__assigned)       
                  ii  = candidate_lines(i).configurations(j).inds_intrsctns__assigned(h);                        
                  intesetcion_coord3D  = intersections(ii).coordinates3D;
                  text(intesetcion_coord3D(1,1), intesetcion_coord3D(1,2), intesetcion_coord3D(1,3), sprintf('%d', ii));
                  plot3( intesetcion_coord3D(1,1), intesetcion_coord3D(1,2), intesetcion_coord3D(1,3),...
                    'g*');
                end 
                for h = 1:length(candidate_lines(i).configurations(j).inds_intrsctns__mult_cnddts)
       
                        ii  = candidate_lines(i).configurations(j).inds_intrsctns__mult_cnddts(h);
                        hi  = candidate_lines(i).configurations(j).inds_intrsctns__mult_cnddts_ind(h);
       
                  intesetcion_coord3D  = intersections(ii).cnddts3D(hi).coordinates3D;
                  text(intesetcion_coord3D(1,1), intesetcion_coord3D(1,2), intesetcion_coord3D(1,3), sprintf('%d_%d', ii, hi));
                  plot3( intesetcion_coord3D(1,1), intesetcion_coord3D(1,2), intesetcion_coord3D(1,3),...
                    'r*');
                end
                disp('');
            end
            
            %Compute coverage:
            [candidate_lines(i).configurations(j).p_coverage,...
             candidate_lines(i).configurations(j).inds_intrsctns,...
             candidate_lines(i).configurations(j).p_intrsctns_dists] = ...
                    computeCost2DCoverage(...
                                candidate_lines(i).configurations(j).inds_intrsctns,...
                                candidate_lines(i).configurations(j).p_intrsctns_dists, ...
                                intersections,...
                                cur_stroke.primitive_geom,...
                                cur_stroke.length2DLikely);%cur_stroke.length2D);

            %Assign directional:
            candidate_lines(i).configurations(j).p_directional = pDirectional{i}(j);

            %Assign directional:
            candidate_lines(i).configurations(j).p_full = ...
                            costComposite(candidate_lines(i).configurations(j).p_coverage, ...
                                          candidate_lines(i).configurations(j).p_directional ,...
                                          cur_stroke.line_group);
             
%             fprintf('Line %d, configuration %d, cost full %.3f \n', i, j, candidate_lines(i).configurations(j).p_full);
               

            if ESTIMATE_JOINTLY    
                configuration = candidate_lines(i).configurations(j);

                cost = configuration.p_full;
                strks_jnt = cur_stroke.ind;
                strks_jnt_cnddt_lns = i;
                strks_jnt_cnddt_lns_cnfgrtns = j;
                
                if (cost > 0.5)
                    try

                    [cost_joint,...
                     strokes_topology,...
                     ~] = ...
                        computeJointCostDependentStrokes(cost,...
                                              strks_jnt,...
                                              strks_jnt_cnddt_lns,...
                                              strks_jnt_cnddt_lns_cnfgrtns,...
                                              configuration, ...
                                              strokes_topology, ...
                                              intersections);
                    catch e
                       rethrow(e); 
                    end
            
                    p_joint{i}(j) = cost_joint;                             
                else
                    p_joint{i}(j) = cost;                             
                end

                candidate_lines(i).configurations(j).p_full_joint = p_joint{i}(j);                        
            end
        end
       [p_joint_cnddt_lns(i), ind_m] = max(p_joint{i}(:));
       p_full_cnddt_lns(i) = candidate_lines(i).configurations(ind_m).p_full;
       candidate_lines(i).max_cost =  p_joint_cnddt_lns(i);
       candidate_lines(i).max_cost_stroke = p_full_cnddt_lns(i);
     
    end
    
    for i = 1:length(hi)
        close(hi{i});
    end
end