function [length_commom, length_s1, length_s2, s_new_points2D,dists] = findEndpointsCommonRegionPolyPoly(s1_points2D, s2_points2D, r)
ind1 = 1:length(s1_points2D);
ind2 = length(s1_points2D)+[1:length(s2_points2D)];

x = [[s1_points2D(:).x] [s2_points2D(:).x]]';
y = [[s1_points2D(:).y] [s2_points2D(:).y]]';

% Find the alrgest diemnsion:
[~,d] = max([max(x) - min(x), max(y) - min(y)]);

% global filepath_sketch_img;
% img = readSketchImg(filepath_sketch_img);

% figure(13); 
% hold off; 
% imshow(img);
% hold on;
% plot([s1_points2D(:).x], [s1_points2D(:).y], 'g:', 'LineWidth', 2);
% plot([s2_points2D(:).x], [s2_points2D(:).y], 'b:', 'LineWidth', 2);

if d == 1
    [y1,x1,ind,dists,er] = fitPoly2(x,y);
    
    if er > r*2
        [x_,y_,ind_,dists_,er_] = fitPoly2(y,x);
        if er_ < er
            y = y_;
            x = x_;
            ind = ind_;
            dists = dists_;
          else
            y = y1;
            x = x1;            
        end
    else
         y = y1;
         x = x1;  
    end
    
else
    [x1,y1,ind,dists,er] = fitPoly2(y,x);   
    if er > r*2
        [y_,x_,ind_,dists_,er_] = fitPoly2(x,y);
        if er_ < er
            y = y_;
            x = x_;
            ind = ind_;
            dists = dists_;
        else
            y = y1;
            x = x1;            
        end
    else
         y = y1;
         x = x1;  
    end
end
% 

% plot(x, y, 'r:', 'LineWidth', 2);

axis equal;

ind1 = find(ismember(ind,ind1));
ind2 = find(ismember(ind,ind2));

length_s1 = lengthStroke( [x(ind1(1):ind1(end)) y(ind1(1):ind1(end))] );
length_s2 = lengthStroke( [x(ind2(1):ind2(end)) y(ind2(1):ind2(end))] );

% dists = dists([ind1(1) ind1(end) ind2(1) ind2(end)]);

v1 = max(ind1(1),ind2(1));
v2 = min(ind1(end),ind2(end));
if v2>v1
    length_commom = lengthStroke([x(v1:v2) y(v1:v2)]);
else
    length_commom =0;
end

num_points = length(x);

s_new_points2D = struct('x', cell(1, num_points), 'y', cell(1, num_points));

xc = num2cell(x);
[s_new_points2D.x] = xc{:};

yc = num2cell(y);
[s_new_points2D.y] = yc{:};


end


function [y,x,ind,dists,er] = fitPoly2(x,y)
[x, ind] = sort(x);
y = y(ind);
% plot(x,y, 'm');
[f, dof] = fit(x, y, 'poly2');
num_points = length(x);
er = dof.rmse;
s_new_points2D = struct('x', cell(1, num_points), 'y', cell(1, num_points));

dists = abs((y-f(x)));
y = f(x);

end