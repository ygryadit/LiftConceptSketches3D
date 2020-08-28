%https://www.intmath.com/applications-differentiation/8-radius-curvature.php
%Input:
%       stroke_points:

function [curvature, G0_points] = polylineCurvature2DFitCircleArrays(stroke_points, indices_sampled)
%     % Compute max dist between points to ensure that always there are at
%     % least three poins.
%     
%     distances = sqrt(sum((stroke_points(1:(end-1),:) - stroke_points(2:end,:)).^2, 2));
    
%     dist_pix = max(15, ceil(max(distances)))
    G0_points = [];
    dist_pix = 6;
    curvature = zeros(1, length(indices_sampled));
    j = 0;
    for p = indices_sampled
        j = j+1;
        try
            Idx = rangesearch(stroke_points,stroke_points(p,:),dist_pix);
        catch
           disp(''); 
        end
%         curvature(p) = computeLocalCurvatureFitCircle(stroke_points(p-1,:), stroke_points(p,:), stroke_points(p+1,:));
%         disp(p)
        
        Idx{1} = min(Idx{1}):max(Idx{1});
        
        if length(Idx{1}) < 4
%            if p == size(stroke_points,1) && min(Idx{1}) ~= 1
%                 Idx{1} = [p-2, p-1, p];
%            elseif p == 1 && max(Idx{1}) ~= size(stroke_points,1)
%                 Idx{1} = [1, 2, 3];    
%            elseif min(Idx{1}) ~= 1 && max(Idx{1})~=size(stroke_points,1)
%                 Idx{1} = [p-1, p, p+1];    
%            end
             curvature(j) = 0.0;
             Par(4) = 0.0;
        else

            [Par, is_singular] = CircleFitByPratt(stroke_points(sort(Idx{1}),:));
           
            curvature(j) = 1.0/Par(3);
            if is_singular
                curvature(j) = 0.0;
            end
        end
        
            
      
%         disp(curvature(p));
%         text(stroke_points(p,1)+4, stroke_points(p,2)+4, num2str(curvature(j)));
        if isnan(curvature(j))
            curvature(j) = 0.0;
        end
%         fprintf('Par(4) = %.2f\n', Par(4));
%         plot(stroke_points(p,1), stroke_points(p,2), 'o');
%         plotCircle(Par(1), Par(2), Par(3), [0,0,0]);   
%         plot(stroke_points( Idx{1},1), stroke_points( Idx{1},2), 'r*');
        if Par(4) > 0.5 && p > (min(Idx{1})) && p < (max(Idx{1}))
           %Check for G_0 discontiniuty
%            plotCircle(Par(1), Par(2), Par(3), [0,0,0]);   
           idx_p = sort(Idx{1});
           i_c = find(idx_p == p);
           
           if length(idx_p(1:i_c)) <= 3
               Par_l(4) = 0;
           
           else
               try
                    Par_l = CircleFitByPratt(stroke_points(idx_p(1:i_c),:));
               catch
                   disp(''); 
               end
           end
           
           if length(idx_p(i_c:end)) <= 3
               Par_r(4) = 0;
           else
               try
               Par_r = CircleFitByPratt(stroke_points(idx_p(i_c:end),:));
               catch
                   disp('');
               end
           end
%             plot(stroke_points(idx_p(1:i_c),1),stroke_points(idx_p(1:i_c),2), 'r*');
%            plotCircle(Par_l(1), Par_l(2), Par_l(3), [1,0,0]);  
%             plot(stroke_points(idx_p(i_c:end),1),stroke_points(idx_p(i_c:end),2), 'g*');
%            plotCircle(Par_r(1), Par_r(2), Par_r(3), [0,1,0]);   
%            if ~isinf(Par_l(3)) && ~isinf(Par_r(3)) && Par_l(4) < 1e-2 && Par_r(4) < 1e-2 && Par(4) >  Par_r(4) &&  Par(4) >  Par_l(4)
%             if Par_l(4) < 1e-2 && Par_r(4) < 1e-2 && Par(4) >  Par_r(4) &&  Par(4) >  Par_l(4)
%                G0_points = [G0_points,j];
%            end
            if Par_l(4) < 0.5 && Par_r(4) < 0.5 
               G0_points = [G0_points,j];
           end
        end    
    end
%     curvature = curvature/num_vals;
curvature = curvature';

end

function res = pointsAreEqual(p1, p2)
    res = p1(1,1) == p2(1,1) || p1(1,2) == p2(1,2);
end

function curvature = computeLocalCurvatureFitCircle(point1, point2, point3)
    m1 = (point2(:,2) - point1(:,2))/(point2(:,1) - point1(:,1));
    m2 = (point3(:,2) - point2(:,2))/(point3(:,1) - point2(:,1));
    
    nom     = m1*m2*(point1(:,2) - point3(:,2)) + m2*(point1(:,1) + point2(:,1)) - m1*(point2(:,1) + point3(:,1));
    denom   = 2*(m2-m1);
    x_c  = nom/denom;
    y_c = -1/m1*(x_c - (point1(:,1) + point2(:,1))/2.0) + (point1(:,2) + point2(:,2))/2.0 ;
    
    r = sqrt((point1(:,1) - x_c)^2+(point1(:,2) - y_c)^2);
    curvature = 1/r;
end


