function [sketch, img, scale] = scaleCenterSketch(sketch, img)
global filepath_sketch_img;
global sketch_height;

new_dim = 512;
sketch_height = new_dim;
% figure(1);
% imshow(img);
% hold on;

    
[min_x, max_x, min_y, max_y] = findBB(sketch);

[center, scale] = findCenterScale(min_x,...
                                  max_x,...
                                  min_y,...
                                  max_y,...
                                  new_dim);

[img, dx, dy] = scaleCenterBitmap(img,...
                                min_x, max_x,...
                                min_y, max_y,... 
                                new_dim, scale);
                            
sketch = scaleCenterVector(sketch, scale,  dx, dy, min_x,  min_y);

sketch.canvas.height = new_dim;
sketch.canvas.width = new_dim;

end

function [min_x, max_x,min_y, max_y] = findBB(sketch)
    max_x = -Inf;
    max_y = -Inf;
    min_x = Inf;
    min_y = Inf;
    for i = 1:length(sketch.strokes)        
        x = [sketch.strokes(i).points(:).x];
        max_x = max([max_x x]);
        
        x = [sketch.strokes(i).points(:).x];
        min_x = min([min_x x]);
        
        y = [sketch.strokes(i).points(:).y];
        max_y = max([max_y y]);
        
        y = [sketch.strokes(i).points(:).y];
        min_y = min([min_y y]);
    end
    
    min_x = (max(floor(min_x-1), 1.0));
    max_x = (min(ceil(max_x+1), sketch.canvas.height));
    min_y = (max(floor(min_y-1), 1.0));
    max_y = (min(ceil(max_y+1), sketch.canvas.height));
%     plot([min_x max_x], [min_y min_y],'r')
%     plot([min_x max_x], [max_y max_y],'r')
%     plot([min_x min_x], [min_y max_y],'r')
%     plot([max_x max_x], [min_y max_y],'r')
end


function [center, scale] = findCenterScale(min_x,...
                                           max_x,...
                                           min_y,...
                                           max_y,...
                                           new_dim)
    
    center = 0.5*[max_x + min_x, max_y + min_y];    
    
    new_max_dim = new_dim - ceil((new_dim - 0.8*new_dim)*0.5)*2;
    
    max_dim = max([(max_x - min_x+1), (max_y - min_y+1)]);
    scale = new_max_dim/max_dim;

end

function sketch = scaleCenterVector(sketch, scale, off_x, off_y, min_x, min_y)
       
    for i = 1:length(sketch.strokes)
        for j = 1:length( sketch.strokes(i).points) 
            sketch.strokes(i).points(j).x = (sketch.strokes(i).points(j).x - min_x)*scale + off_y;
            sketch.strokes(i).points(j).y = (sketch.strokes(i).points(j).y - min_y)*scale + off_x;
           
        end
%         plot([sketch.strokes(i).points(:).x],[[sketch.strokes(i).points(:).y]]);
    end

end

function [img_new, off_x, off_y] = scaleCenterBitmap(img,...
                                min_x, max_x,...
                                min_y, max_y,... 
                                new_dim, scale)

   img = img(min_y:max_y,min_x:max_x,:);
   img = imresize(img, scale);
%    figure(2); imshow(img);
   img_new = 255*ones(new_dim, new_dim, 3, 'uint8');
   
   off_x_ = 0.5*(new_dim - size(img,1));   
   off_x = floor(off_x_);
   off_y_ = 0.5*(new_dim - size(img,2));
   off_y = floor(off_y_);
   
   dx = off_x_ - off_x;
   dy = off_y_ - off_y;
   
   img_new(off_x:(off_x+size(img,1)-1),off_y:(off_y+size(img,2)-1),:) = img; 
%    figure(3); imshow(img_new);
%    hold on;
end
