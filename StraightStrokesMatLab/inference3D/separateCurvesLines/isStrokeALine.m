function [is_line, line] = isStrokeALine(x,y)
    is_line = false;
    
%     [x, ind_unique, ~] = unique(x, 'stable');
%     y = y(ind_unique);
%     [y, ind_unique, ~] = unique(y, 'stable');
%     x = x(ind_unique);
    
    if length(x) < 3
        line = [];
        return;
    end
    
%     mean_x = mean(x);
%     mean_y = mean(y);
%     
%     plot(mean_x, mean_y,'r*', 'LineWidth', 2);
    

    % Compute weigthed mean:
        x_ = [x(1); x; x(end)];
        y_ = [y(1); y; y(end)];
    
        seg_lengths = sqrt((y_(1:end-1)- y_(2:end)).^2 + (x_(1:end-1)- x_(2:end)).^2);
        weights = seg_lengths(1:end-1) + seg_lengths(2:end);
        mean_x = sum(weights.*x)/sum(weights);
        mean_y = sum(weights.*y)/sum(weights);
    
%         plot(mean_x, mean_y,'b*', 'LineWidth', 2);
    %
    
    zmx = (x-mean_x);
    zmy = (y-mean_y); 
%     D = [sum(zmx.^2) sum(zmx.*zmy); sum(zmx.*zmy) sum(zmy.^2)];
    D = [sum(zmx.^2.*weights) sum(zmx.*zmy.*weights); sum(zmx.*zmy.*weights) sum(zmy.^2.*weights)];
    
    [~, L] = eig(D);

    
%     if (L(1,1)/L(2,2) > 5e-4)
    if (L(1,1)/L(2,2) > 1e-3)
        line = [];
        return;
    end
    
    is_line = true;
    

%     
%    A = dir ;
%    b = [x(1) - mean_x, y(1) - mean_y];
%    t1 = A/b;
%    
%    x1 = mean_x + t1*dir(1);
%    y1 = mean_y + t1*dir(2);     
%    
%    
%    b = [x(end) - mean_x, y(end) - mean_y];
%    t2 = A/b;
%    
%    x2 = mean_x + t2*dir(1);
%    y2 = mean_y + t2*dir(2);     
   
    % Find longest dimension:
    [~, ind] = max([max(x)-min(x) max(y)-min(y)]);
    if (ind == 1)
        [f, ~] = fit(x,y,'poly1');
    else
        [f, ~] = fit(y,x,'poly1');
    end
%     error = gof.rmse;

    % Normilise to unit area
%     scale = (max(x)-min(x))*(max(y)-min(y));

    if (ind == 1)
        [x1, ind_x1] = min(x);
        [x2, ind_x2] = max(x);
        y1 = f(x1);
        y2 = f(x2);  
        
        % Ensure that the line endpoints are in the same order as in the
        % polyline:
        if (ind_x1 < ind_x2)
            line = [x1 x2 y1 y2];
        else
            line = [x2 x1 y2 y1];
        end
    else
        [y1, ind_y1] = min(y);
        [y2, ind_y2] = max(y);
        x1 = f(y1);
        x2 = f(y2);
        
        % Ensure that the line endpoints are in the same order as in the
        % polyline:
        if (ind_y1 < ind_y2)
            line = [x1 x2 y1 y2];
        else
            line = [x2 x1 y2 y1];
        end
    end

    
    
%    dir = [V(1,2), V(2,2)];
%    
%    v = [x1 - mean_x, y1 - mean_y];
%    
%    p1 = [mean_x, mean_y] + (v*dir')/(dir*dir')*dir;
%    
%    v = [x2 - mean_x, y2 - mean_y];
%    
%    p2 = [mean_x, mean_y] + (v*dir')/(dir*dir')*dir;
%     
%     
%     x1 = p1(1);
%     x2 = p2(1);
%     y1 = p1(2);
%     y2 = p2(2);

    
    
    
    
%     plot([x1 x2], [y1 y2], ':');
    
       
% %     length = sqrt((x2-x2)^2 + (f(x2)-f(x1))^2);
%     length = sqrt((x2-x2)^2 + (y2-y1)^2);
% %     error = error/scale;
%     error = error/length;
%     
% %     if error < 1e-5
%     if error < 0.05
%         is_line = true;
%     end
%     
  
    
%     theta = atan2(V(2, 2) , V(1, 2));                       
%    
%     r = sqrt((max(x)-min(x))^2 + (max(y)-min(y))^2);
%     x1 = mean_x - cos(theta)*r/2;
%     x2 = mean_x + cos(theta)*r/2;
%     y1 = mean_y - sin(theta)*r/2;
%     y2 = mean_y + sin(theta)*r/2;            

%     r = mean_x*cos(theta)+mean_y*sin(theta);
    
%     if isnan(y1) || isnan(y2)
%         disp('check me');
%     end

%     theta = NaN;
%     rho = NaN;
%     line = [x1 x2 y1 y2 theta rho];
  
 
%     figure(2);
%     if (ind == 1)
%         plot(f,x,y);
%         hold on;
%         plot(line(1:2), line(3:4), ':');
%         hold off;
%     else
%         plot(f,y,x);    
%         hold on;
%         plot(line(1:2), line(3:4), ':');
%         hold off;
%     end     
%     axis equal;
%     axis ij;
%     
%     x = x - mean(x);
%     y = y - mean(y);
%     D = [sum(x.^2) dot(x,y); dot(x,y) sum(y.^2)];
%     lambda = eig(D);
%     
%     legend('data', sprintf('er = %.5f, r = %.5f', error, lambda(1)/lambda(2)));
end