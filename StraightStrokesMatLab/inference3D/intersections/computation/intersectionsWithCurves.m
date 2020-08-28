function [intersections] = intersectionsWithCurves(strokes_topology, img)
    intersections = [];
    ii = 0;
    debug_cci = false;
    % I. Extend all curves by 2*accuracy radius and strokes by accuracy radius
    
%     strokes_topology = extendStrokes(strokes_topology, img);
    
    % II. Compute curve intersections with straigth strokes
    inds_curves = find(cat(1, strokes_topology(:).primitive_type) == 1);
    inds_lines  = find(cat(1, strokes_topology(:).primitive_type) == 0);
    inds_curves = reshape(inds_curves, 1, []);
    inds_lines = reshape(inds_lines, 1, []);
    
    for ic = inds_curves
        for il  = [inds_lines, setdiff(inds_curves, ic)]
            xv1 = strokes_topology(ic).poly2d_extended(:,1)';
            yv1 = strokes_topology(ic).poly2d_extended(:,2)';
            xv2 = strokes_topology(il).poly2d_extended(:,1)';
            yv2 = strokes_topology(il).poly2d_extended(:,2)';
            
            try
                [d_min, xc_is, yc_is, idxc_is, is_vertex] = poly_poly_dist(xv1, yv1, xv2, yv2);
            catch
                continue;
            end
        
            ind_s2 = max(ic,il);
            idxc_is(idxc_is < 1) = idxc_is(idxc_is < 1) + 1;
            if d_min < strokes_topology(ind_s2).accuracy_radius
                    
                    for k = 1:size(xc_is,1)
                        ii = ii + 1;
                        intersections.coordinates2D(ii,:) = [mean(xc_is(k,:)) mean(yc_is(k,:))];
                        
                        
                        [t1, p1] = findPointPorjectedPositionOnTheSegment([xv1(idxc_is(k,1)),xv1(idxc_is(k,1)+1),...
                                                                  yv1(idxc_is(k,1)),yv1(idxc_is(k,1)+1)],...
                                                                  intersections.coordinates2D(ii,:));
                    
                        t1 = t1 + idxc_is(k,1); 
                    
                        [t2, p2] = findPointPorjectedPositionOnTheSegment([xv2(idxc_is(k,2)),xv2(idxc_is(k,2)+1),...
                                                                       yv2(idxc_is(k,2)),yv2(idxc_is(k,2)+1)],...
                                                                       [mean(xc_is),mean(yc_is)]);
                        t2 = t2 + idxc_is(k,2); 
                     
                        
                        dir1 = strokes_topology(ic).poly2d_extended(idxc_is(k,1)+1,:) - ....
                                strokes_topology(ic).poly2d_extended(idxc_is(k,1),:);
                        dir2 = strokes_topology(il).poly2d_extended(idxc_is(k,2)+1,:) - ....
                                strokes_topology(il).poly2d_extended(idxc_is(k,2),:);
                   
                        cos_dirs = dot(dir1./norm(dir1), dir2./norm(dir2));
                        if strokes_topology(il).primitive_type == 0
                            intersections.tangent(ii) = abs(cos_dirs) > cosd(10);
                            intersections.collinear(ii) = false;
                        else
                            
                            intersections.collinear(ii) = abs(cos_dirs) > cosd(35);
                            intersections.tangent(ii) = false;
                            
%                             plot(strokes_topology(ic).poly2d_extended([idxc_is(k,1) idxc_is(k,1)+1],1), ...
%                                  strokes_topology(ic).poly2d_extended([idxc_is(k,1) idxc_is(k,1)+1],2),'*r');
%                             
%                             plot(strokes_topology(il).poly2d_extended([idxc_is(k,2) idxc_is(k,2)+1],1), ...
%                                  strokes_topology(il).poly2d_extended([idxc_is(k,2) idxc_is(k,2)+1],2),'*g');
%                             
%                              if intersections.collinear(ii)
%                                  plot(intersections.coordinates2D(ii,1),intersections.coordinates2D(ii,2),'o');
%                              end
%                             
%                              disp('');
                        end
                        intersections.strokes_indices(ii,:) = [ic, il];
                        intersections.seg_nums(ii,:) = [t1, t2];
                        intersections.p_dist_str_segs(ii,:) = 1.0;
                        
                        
                        if debug_cci
                            figure(15);
                            hold off;
                            imshow(img);
                            hold on;
                            plot(xv1,yv1,'r');
                            plot(xv2,yv2,'g');

                            plot(xc_is,yc_is,'*c');


                            plot(xv1(idxc_is(:,1)),yv1(idxc_is(:,1)),'*r');
                            plot(xv1(idxc_is(:,1)+1),yv1(idxc_is(:,1)+1),'*r');
                            plot(p1(1),p1(2),'or');


                             plot(xv2(idxc_is(:,2)),yv2(idxc_is(:,2)),'*r');
                             plot(xv2(idxc_is(:,2)+1),yv2(idxc_is(:,2)+1),'*r');
                             plot(p2(1),p2(2),'og');
                         end  
                        
                        
                    end
            end
        end
        
    end
    
    % III. Mark intersections that are likely to be a tangential constraint
    % This should nto be accounted for than computing likely with strokes.

end



