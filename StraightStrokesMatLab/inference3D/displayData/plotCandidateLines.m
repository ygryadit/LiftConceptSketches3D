function plotCandidateLines(candidate_lines, ...
                            strokes_topology, ...
                            cur_stroke,...
                            intersections, fig_num, plot_intersections)
 
    global ESTIMATE_JOINTLY;
    if ~exist('plot_intersections', 'var')
        plot_intersections = true;
    end
    figure(fig_num);
    close(fig_num);
       plotStrokesTopology(fig_num, strokes_topology( cat(1, strokes_topology(:).depth_assigned)));
%     plotStrokesTopolgyIntersectionsTypes(...
%                             strokes_topology( cat(1, strokes_topology(:).depth_assigned)),...
%                             cur_stroke,...
%                             intersections,...
%                             fig_num);
      

    for i = 1:length(candidate_lines)
%        if abs(sum( candidate_lines(i).coordinates3D(1,1:3) - candidate_lines(i).coordinates3D(1,4:6))) > 1e-2
%             plot3(  candidate_lines(i).coordinates3D(1,[1,4]),...
%                 candidate_lines(i).coordinates3D(1,[2,5]),...
%                 candidate_lines(i).coordinates3D(1,[3,6]),...
%                 ':', 'LineWidth', 2);
%        else
           plot3(   candidate_lines(i).coordinates3D_prior(1,[1,4]),...
                    candidate_lines(i).coordinates3D_prior(1,[2,5]),...
                    candidate_lines(i).coordinates3D_prior(1,[3,6]),...
                    ':', 'LineWidth', 2);
                
            if ~plot_intersections
                continue;
            end
            
           for j = 1:length( candidate_lines(i).configurations) 
              
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
                
                  % Find a relevant candidate line and plot it:
                  intrsctng_stk = setdiff(intersections(ii).strokes_indices, cur_stroke.ind);
                  mhs = find(intersections(ii).strokes_indices == intrsctng_stk);
                  additionl_lines  = intersections(ii).cnddts3D(hi).cnddt_lns{mhs};
                  for ki = 1:length(additionl_lines)
                      lmc = additionl_lines(ki);
                      plot3( strokes_topology(intrsctng_stk).candidate_lines(lmc).coordinates3D_prior(1,[1,4]),...
                             strokes_topology(intrsctng_stk).candidate_lines(lmc).coordinates3D_prior(1,[2,5]),...
                             strokes_topology(intrsctng_stk).candidate_lines(lmc).coordinates3D_prior(1,[3,6]),...
                        '.-', 'LineWidth', 2);
                      
                        inds_confgrtns = intersections(ii).cnddts3D(hi).cnfgrtns{mhs}{ki};
                      
                        %intersections assigned strokes:
                        try
                        intrsctns_assgnd = ...
                            [strokes_topology(intrsctng_stk).candidate_lines(lmc).configurations(inds_confgrtns).inds_intrsctns__assigned];
                        catch e
                            rethrow(e)
                        end
                        if ~isempty(intrsctns_assgnd)
                            for lia = intrsctns_assgnd
                                  intesetcion_coord3D  = intersections(lia).coordinates3D;
                                  text(intesetcion_coord3D(1,1), intesetcion_coord3D(1,2), intesetcion_coord3D(1,3), sprintf('%d', lia));
                                  plot3( intesetcion_coord3D(1,1), intesetcion_coord3D(1,2), intesetcion_coord3D(1,3),...
                                    'b*');
                            end
                        end
                        
                        % intersections strokes multiple candidates:
                        intrsctns_mcndtts = cat(2, ...
                            strokes_topology(intrsctng_stk).candidate_lines(lmc).configurations(inds_confgrtns).inds_intrsctns__mult_cnddts);
                        intrsctns_mcndtts_vrsns = cat(2, ...
                            strokes_topology(intrsctng_stk).candidate_lines(lmc).configurations(inds_confgrtns).inds_intrsctns__mult_cnddts_ind);
                        
                        for ia = 1:length(intrsctns_mcndtts)
                              lia = intrsctns_mcndtts(ia);
                              hi = intrsctns_mcndtts_vrsns(ia);
                              try
                                intesetcion_coord3D  = intersections(lia).cnddts3D(hi).coordinates3D;
                              catch
                                dips('');  
                              end
                              text(intesetcion_coord3D(1,1), intesetcion_coord3D(1,2), intesetcion_coord3D(1,3), sprintf('%d_%d', lia, hi));
                              plot3( intesetcion_coord3D(1,1), intesetcion_coord3D(1,2), intesetcion_coord3D(1,3),...
                                'b*');
                        end
                        
                        
                  end
                  
               end
           end
%        end
    
    end
    
%     
%     if ESTIMATE_JOINTLY
%        num_strokes_dependent = sum(~cur_stroke.mask_indcs_intrsctns_prvs_strks_actv);
% 
%        inter_mult_hypothesis = cur_stroke.inds_intrsctns_prvs_strks(~cur_stroke.mask_indcs_intrsctns_prvs_strks_actv);
%        strokes_mult_hypothesis = cur_stroke.indcs_intrsctng_prvs_strks(~cur_stroke.mask_indcs_intrsctns_prvs_strks_actv);
% 
%        for i = 1:num_strokes_dependent
%            i_stroke = strokes_mult_hypothesis(i);
%            if ~isfield(strokes_topology(i_stroke), 'candidate_lines')
%                 continue;
%            end
% %                    fprintf('Number of hypothesis %d\n', length(strokes_topology(i_stroke).candidate_lines));
%            for hj = 1:strokes_topology(i_stroke).num_candidate_lines
%  
%                intesetcion_coord3D  = intersections(inter_mult_hypothesis(i)).candidates3D(hj).coordinates;
%                text(intesetcion_coord3D(1,1), intesetcion_coord3D(1,2), intesetcion_coord3D(1,3),num2str(inter_mult_hypothesis(i)));
%                plot3( intesetcion_coord3D(1,1), intesetcion_coord3D(1,2), intesetcion_coord3D(1,3),...
%                     'g*');
%            end
%        end
%    end
  hf = figure(fig_num);
  set(hf, 'Position',  [749.0000  914.6000  846.4000  647.2000]);
    set(gca, 'CameraPosition', [-1.6560 -3.3123 -1.6620]);
    set(gca, 'CameraTarget', [-0.1179 0.0548 0.1537]);
    
%     saveMP4s(strokes_topology);
end

function plotStrokesTopology(fig_num, strokes_assigned) 
%     figure(fig_num);
%     holf off;
    num_strokes = length(strokes_assigned);
    
    colors = uint8(colormap(parula(num_strokes))*255);
   
    h = figure(fig_num);
     
%     set(h,'Position', [1168 74 1386 1195]);
%     set(h,'Position', [3522         178         958        1074]);
%     hold off;
%    
   
%         view(-30.2000, 21.2000);
    set(h, 'Name', '3D candidate intersections');
   

    
    for i = 1:length(strokes_assigned)
        plot3(  strokes_assigned(i).points3D(:,1),...
                strokes_assigned(i).points3D(:,2),...
                strokes_assigned(i).points3D(:,3), ...
                'LineWidth', 1.0,...
                'Color', colors(i, :));
         hold on;
    end
    grid on;
    hold on;
    h = colorbar;
     
    b = num2str([1:num_strokes]');
    c = cellstr(b);
    h.Ticks = [0:(num_strokes-1)]/(num_strokes-1);
    h.TickLabels = c;
    
     axis equal;
    global ZUP
    if ~ZUP
        axis ij
        set(gca, 'Zdir', 'reverse')
    end
    global CameraView;
    view(3);
%     view(CameraView(1),CameraView(2));
%     camva(11.2020);
%     campos([-1.0327   -0.7573   -0.5696]);
%     camtarget([0.0319   0.0426    0.1172]);
%     camup([0     0    -1]);
end