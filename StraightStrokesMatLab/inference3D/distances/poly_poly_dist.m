function [d_min, xc_is, yc_is, idxc_is, is_vertex] = poly_poly_dist(xv1, yv1, xv2, yv2)
%poly_poly_dist Find minimum distance between two polylines.
%
% Description:
% Polyline is defined as a set of nv-1 segments connecting nv ordered 
% vertices v(1), v(2), ..., v(nv).
% Distance between 2 non-intersecting polylines poly1 and poly2 is defined
% as a minimum of distances of poly1 vertices to poly2 segments, or poly2
% segments to poly1 vertices.
% If polylines intersect, the distance is set to zero, and the coordinates
% of intersection points are returned.
%
% Input arguments:
% xv1 - vector of X coordinates of the first polyline vertices (1 X nv1)
% yv1 - vector of Y coordinates of the first polyline vertices (1 X nv1)
% xv2 - vector of X coordinates of the second polyline vertices (1 X nv2)
% yv2 - vector of Y coordinates of the second polyline vertices (1 X nv2)
%
%
% Output arguments:
%
% d_min - distance between 2 polyline (scalar). If polylines interesect,
% d_min is set to zero
% 
% xc_is - X coordinates of closest points (nc_is X 2), where nc_is is a 
% number of closest or intersection points. nc_is can be greater than 1 in 
% case when there are more than 1 points with SAME minimum distance. 
% xc_is(:,1)contains X coordinates of closest or intersection points 
% belonging to the first polyline, xc_is(:,2) - to the second polyline. If 
% polylines intersect (d_min==0), xc_is contains the coordinates of 
% intersection points, and xc_is(:,1) = xc_is(:,2).
%
% yc_is - Y coordinates of closest or intersection points (nc_is X 2)
%
% idxc_is - indices of polyline segments that contain the closest or 
% interesection points (nc_is X 2). idxc_is(:,1) contains the indices of 
% first polyline's segments, idxc_is(:,2) - of second polyline's segments.
%
% is_vertex - logical array (nc_is X 2). If is_vertex(j,1) is true, the
% closest or intersection point j that belongs to poly1 is a vertex. If 
% is_vertex(j,2) is true, the closest or intersection point j that belongs
% to poly2 is a vertex
%
% Revision history:
% Oct 28, 2015 - Created (Michael Yoshpe).
%**************************************************************************
% some inputs sanity checks
nv1 = length(xv1);
nv2 = length(xv2);
if(nv1 < 2)
   error('Polyline 1 must have at least 2 vertices');
end
if(nv2 < 2)
   error('Polyline 2 must have at least 2 vertices');
end
% number of vertices in each polyline
xmin1 = min(xv1);
xmax1 = max(xv1);
ymin1 = min(yv1);
ymax1 = max(yv1);
xmin2 = min(xv2);
xmax2 = max(xv2);
ymin2 = min(yv2);
ymax2 = max(yv2);
% distance from vertices of poly1 to segments or vertices of poly 2
[d_min12, x_d_min12, y_d_min12, is_vertex12, idx_c12, xc12, yc12, is_in_seg12, Cer12, Ppr12] = ...
   p_poly_dist(xv1, yv1, xv2, yv2, false);
% crude test for intersection of polylines
if((xmax1 < xmin2) | (xmin1 > xmax2) | (ymax1 < ymin2) | (ymin1 > ymax2))
   % polylines definitely don't intersect, set intersection flag to false
   is_flag = false;
else
   % polyline could still intersect, check for possible intersection
   % points of poly1 in rotated coordinate systems defined by segments of 
   % poly2
   xpr = zeros(nv1, (nv2-1));
   ypr = zeros(nv1, (nv2-1));
   xpr(:,:) = Ppr12(1,:,:);
   ypr(:,:) = Ppr12(2,:,:);
   idx_is = (sign(ypr(1:(end-1),:)) ~= sign(ypr(2:end,:)));
   
   if(any(any(idx_is))) 
      % some vertices of poly1 have opposite Y coordinates, when 
      % transformed to coordinate systems defined by segments of poly2 
      % Continue checking for intersection
      [ii1, jj1] = find(idx_is);
      idxl1 = sub2ind(size(ypr), ii1, jj1);
      idxl2 = sub2ind(size(ypr), ii1+1, jj1);
      
      % x coordinates of a lines connecting poly1 vertices in rotated
      % system. To intersect, they should fall WITHIN poly2 segments
      xr_is = (xpr(idxl1).*ypr(idxl2) - xpr(idxl2).*ypr(idxl1))./(ypr(idxl2) - ypr(idxl1));
      
      % poly2 segments lengths
      vds_is2 = hypot((xv2(jj1+1)-xv2(jj1)), (yv2(jj1+1)-yv2(jj1)));
      ii_is = find((xr_is >= 0) & (xr_is(:) <= vds_is2(:)));
      
      if(~isempty(ii_is)) % intersection points found
         is_flag = true;
         n_is = length(ii_is); % number of segments intersections
         % Polylines definitely interesect, since some segments of
         % poly1 cross segments of poly2. Translate the intersection 
         % points into original coordinate system
         
         Pcr = zeros([2 n_is]);
         Pcr(1, :) = xr_is(ii_is);
         Pcr(2, :) = 0;
         
         % Cre is a rotation matrix from rotated to original system
         Cre = permute(Cer12(:,:,jj1(ii_is)), [2 1 3]);
         
         Cre1 = zeros(2, n_is);
         Cre2 = zeros(2, n_is);
         
         Cre1(:,:) = Cre(1, :, :);
         Cre2(:,:) = Cre(2, :, :);
                  
         % coordinates of consecutive vertices
         P1 = [xv2(:) yv2(:)];
         % Pce is a 2D array of size 2 * n_is that holds the coordinates of
         % intersection points in original coordinate systems. Pce(1,j) is 
         % an X coordinate of the intersection of point, Pce(2,j) is its 
         % Y coordinate
         Pce = zeros(2, n_is);
         Pce(1,:) = Cre1(1,:).*Pcr(1,:) + Cre1(2,:).*Pcr(2,:);
         Pce(2,:) = Cre2(1,:).*Pcr(1,:) + Cre2(2,:).*Pcr(2,:);
                          
         % Adding the P1 vector
         Pce(1,:) =  Pce(1,:) + P1(jj1(ii_is),1)';
         Pce(2,:) =  Pce(2,:) + P1(jj1(ii_is),2)';
         
         % x and y coordinates of the projected (cross-over) points in 
         % original coordinate frame
         d_min = 0;
         xc_is = zeros(n_is, 2);
         yc_is = zeros(n_is, 2);
         xc_is(:,1) = Pce(1,:);
         xc_is(:,2) = Pce(1,:);
         yc_is(:,1) = Pce(2,:);
         yc_is(:,2) = Pce(2,:);
         idxc_is = [ii1(ii_is) jj1(ii_is)];
         
         % If polylines intersect, intersection points are most often lie
         % INSIDE segments, initialize is_vertex to false
         is_vertex = false(n_is, 2);
         
         % check for corner cases, where polylines intersect at vertices
                
         % check if interesection points are close to poly1 vertices
         cond_poly1 = ((abs(xv1(idxc_is(:,1))'-xc_is(:,1)) < 5*eps) & ...
                      (abs(yv1(idxc_is(:,1))'-yc_is(:,1)) < 5*eps)) | ...
                      ((abs(xv1(idxc_is(:,1)+1)'-xc_is(:,1)) < 5*eps) & ...
                      (abs(yv1(idxc_is(:,1)+1)'-yc_is(:,1)) < 5*eps));
         % check if interesection points are close to poly2 vertices                   
         cond_poly2 = (Pcr(1, :) < 10*eps) | (Pcr(1, :) > vds_is2(ii_is)-10*eps);
         
         is_vertex(:, 1) = cond_poly1;
         is_vertex(:, 2) = cond_poly2;
      else
         % polylines don't interesect
         is_flag = false;
      end
   else 
      % polylines don't interesect
      is_flag = false;
   end
   
end
% Plylines don't intersect, find the minimum distance
if(~is_flag)
   % distance from vertices of poly1 to segments or vertices of poly2
   [d_min21, x_d_min21, y_d_min21, is_vertex21, idx_c21, xc21, yc21, is_in_seg21, Cer21, Ppr21] = ...
      p_poly_dist(xv2, yv2, xv1, yv1, false);
   
   % find minimum distance
   [dm12, im12] = min(d_min12);
   [dm21, im21] = min(d_min21);
   
   [d_min, imj] = min([dm12, dm21]);
   
   % find closest points (could be more than 1 for each polyline, with same
   % distance) 
   im1 = find(abs(d_min12-d_min)<5*eps);
   im2 = find(abs(d_min21-d_min)<5*eps);
   
   xc_is1 = [x_d_min21(im2) xv2(im2)';...
            xv1(im1)' x_d_min12(im1)];
         
   yc_is1 = [y_d_min21(im2) yv2(im2)';...
            yv1(im1)' y_d_min12(im1)];
    
   % note that im1 and im2 are indices of closest POINTS. To find the 
   % indices of closest SEGMENTS, we have to substract 1
   im1s = im1;
   im2s = im2;
   
   if(~isempty(im1))
      im1s = im1 - 1;
   end
   if(~isempty(im2))
      im2s = im2 - 1;
   end
   
   idxc_is1 = [idx_c21(im2) im2s;...
              im1s idx_c12(im1)];
   
   % remove possible duplicates        
   [idxc_is, ia, ic] = unique(idxc_is1, 'rows');
   
   xc_is = xc_is1(ia,:);
   yc_is = yc_is1(ia,:);   
   
   is_vertex1 = [is_vertex21(im2) true(length(im2),1);
                true(length(im1),1) is_vertex12(im1)];
    
   is_vertex = is_vertex1(ia, :);          
             
end
end