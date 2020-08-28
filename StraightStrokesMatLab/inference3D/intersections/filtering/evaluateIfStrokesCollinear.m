function collinear = evaluateIfStrokesCollinear(strokes_topology, ind_s1, ind_s2)
    collinear = false;
    xv1 = strokes_topology(ind_s1).poly2d_extended(:,1)';
    yv1 = strokes_topology(ind_s1).poly2d_extended(:,2)';
    xv2 = strokes_topology(ind_s2).poly2d_extended(:,1)';
    yv2 = strokes_topology(ind_s2).poly2d_extended(:,2)';
        
    
   
    [d_min, xc_is, yc_is, idxc_is, is_vertex] = poly_poly_dist(xv1, yv1, xv2, yv2);
   
    idxc_is(idxc_is == 0) = 1;
                  
    for k = 1:size(xc_is,1)


        dir1 = strokes_topology(ind_s1).poly2d_extended(idxc_is(k,1)+1,:) - ....
                strokes_topology(ind_s1).poly2d_extended(idxc_is(k,1),:);
        dir2 = strokes_topology(ind_s2).poly2d_extended(idxc_is(k,2)+1,:) - ....
                strokes_topology(ind_s2).poly2d_extended(idxc_is(k,2),:);

        cos_dirs = dot(dir1./norm(dir1), dir2./norm(dir2));

        if (strokes_topology(ind_s2).primitive_type == 1 & strokes_topology(ind_s1).primitive_type == 0) | ...
            (strokes_topology(ind_s2).primitive_type == 0 & strokes_topology(ind_s1).primitive_type == 1)
            collinear = false;
        elseif strokes_topology(ind_s2).primitive_type == 1 | strokes_topology(ind_s1).primitive_type == 1            

            collinear = collinear | (abs(cos_dirs) > cosd(35));
            
            

    %                             plot(strokes_topology(ind_s1).poly2d_extended([idxc_is(k,1) idxc_is(k,1)+1],1), ...
    %                                  strokes_topology(ind_s1).poly2d_extended([idxc_is(k,1) idxc_is(k,1)+1],2),'*r');
    %                             
    %                             plot(strokes_topology(ind_s2).poly2d_extended([idxc_is(k,2) idxc_is(k,2)+1],1), ...
    %                                  strokes_topology(ind_s2).poly2d_extended([idxc_is(k,2) idxc_is(k,2)+1],2),'*g');
    %                             
    %                              if intersections.collinear(ii)
    %                                  plot(intersections.coordinates2D(ii,1),intersections.coordinates2D(ii,2),'o');
    %                              end
    %                             
    %                              disp('');
        elseif strokes_topology(ind_s2).primitive_type == 0 & strokes_topology(ind_s1).primitive_type == 0            
            collinear = collinear | (abs(cos_dirs) > cosd(5));
        end
       
    end
end