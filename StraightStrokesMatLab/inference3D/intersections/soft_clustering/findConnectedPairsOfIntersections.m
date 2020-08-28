% This function computed pairs of intersections which are within accuracy
% radiuses of each other, meaning that they can be visually grouped.
% 
% 
function [pairsInterInter] = findConnectedPairsOfIntersections(strokes_topology, intersections, img)
    global fid;
    tic
   
    
    tic
    pairsInterInter = findPairsInterInter(intersections, strokes_topology, img);
    elapsed_time = toc;
    fprintf(fid, 'Time to find linked pairs of intersections: %.3f\n', elapsed_time);
    
    pairsInterInter = orderPairs(pairsInterInter);
    
%     tic
%     pairsInterAttrac = findPairsInterAttraction(intersections, points_attraction_to_pair);
%     elapsed_time = toc;
%     fprintf('Time to find linked pairs of attraction points: %.3f\n', elapsed_time);
%     
%      plotPairs(intersections, pairsInterInter, img);
end


function plotPairs(intersections, pairsInterInter, img)
    
    
    for i = 1:length(pairsInterInter)
       figure(21); 
       hold off;
       imshow(img);
       hold on;
       plot(intersections.coordinates2D(pairsInterInter(i,:),1),...
           intersections.coordinates2D(pairsInterInter(i,:),2), '*-');
    end

end

function plotAccuracyRadiusIntersections(intersections, img)
    figure; imshow(img);
    hold on;
    
    for i =1:size(intersections.coordinates2D,1)
       plot(intersections.coordinates2D(i,1),...
           intersections.coordinates2D(i,2), '*');
       plotCircle(intersections.coordinates2D(i,1),...
           intersections.coordinates2D(i,2), intersections.accuracy_radius(i), [0,0,0]);
    end

end


function pairsInterAttrac = findPairsInterAttraction(intersections, points_attraction_to_pair)
      num_intersections = size(intersections.coordinates2D,1);
      num_attractions = size(points_attraction_to_pair.coordinates2D,1);
      [INT, ATTR] = meshgrid(1:num_intersections,1:num_attractions);
      INT = INT(:);
      ATTR = ATTR(:);

      
      distances1 = sqrt(sum((intersections.coordinates2D(INT,:) - points_attraction_to_pair.coordinates2D(ATTR,[1,2])).^2,2));
%       distances2 = sqrt(sum((intersections.coordinates2D(INT,:) - points_attraction_to_pair.coordinates2D(ATTR,[3,4])).^2,2));
       
%       mask_to_pair = (distances1 < intersections.accuracy_radius(INT)) | ...
%                      (distances2 < intersections.accuracy_radius(INT)) | ...
%                      (distances1 < points_attraction_to_pair.accuracy_radius(ATTR)) | ...
%                      (distances2 < points_attraction_to_pair.accuracy_radius(ATTR));
      mask_to_pair = (distances1 < intersections.accuracy_radius(INT)) | ...
                     (distances1 < points_attraction_to_pair.accuracy_radius(ATTR));
                 
%      ind_to_pair = find(mask_to_pair);

     pairsInterAttrac = [INT(mask_to_pair) ATTR(mask_to_pair)];
end
 
 

