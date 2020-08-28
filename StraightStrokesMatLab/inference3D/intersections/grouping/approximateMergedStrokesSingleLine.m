function s_new_points2D = approximateMergedStrokesSingleLine(s1_points2D, s2_points2D)

x = [[s1_points2D(:).x] [s2_points2D(:).x]]';
y = [[s1_points2D(:).y] [s2_points2D(:).y]]';

% Find the alrgest diemnsion:
[~,d] = max([max(x) - min(x), max(y) - min(y)]);

%     global filepath_sketch_img;
%     img = readSketchImg(filepath_sketch_img);
    
if d == 1
    [y,x,f] = fitPoly2(x,y);
    
% 
%     figure(13); 
%     hold off; 
%     imshow(img);
%     hold on;
%     plot(f, x, y);
%     axis equal;
% 

else
    [x,y,f] = fitPoly2(y,x);
    

%     figure(13); 
%     hold off; 
%     imshow(img);
%     hold on;
%     plot(f, y, x);
%     axis equal;

end

num_points = length(x);

s_new_points2D = struct('x', cell(1, num_points), 'y', cell(1, num_points));

xc = num2cell(x);
[s_new_points2D.x] = xc{:};

yc = num2cell(y);
[s_new_points2D.y] = yc{:};


end

function [y,x,f] = fitPoly2(x,y)
[x, ind] = sort(x);
y = y(ind);

f = fit(x, y, 'poly2');
num_points = length(x);

s_new_points2D = struct('x', cell(1, num_points), 'y', cell(1, num_points));

y = f(x);

end