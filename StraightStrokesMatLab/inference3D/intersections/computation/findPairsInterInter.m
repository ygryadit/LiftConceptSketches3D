function [pairsInterInter, distances] = findPairsInterInter(intersections, strokes_topology, img)
    
    pairsInterInter =[];
    for i = 1:length(strokes_topology)
        if strokes_topology(i).primitive_type ~= 0 && ...
           strokes_topology(i).primitive_type ~= 1       
            continue;
        end
       
        %Copmute grouping only for straight strokes and curves.
        %Copmute grouping only along each stroke (simplified version)
        
        [INT1, INT2] = meshgrid(strokes_topology(i).indcs_intrsctns);
        num_intersections = length(strokes_topology(i).indcs_intrsctns);
        A = true(num_intersections);
        A = triu(A,1);
        INT1 = INT1(A);
        INT2 = INT2(A);
        
        distances = sqrt(sum((intersections.coordinates2D(INT1,:) - intersections.coordinates2D(INT2,:)).^2,2));
        
        mask_to_pair =  (distances < intersections.accuracy_radius(INT1)) | ...
                        (distances < intersections.accuracy_radius(INT2));
    
        distances = distances(mask_to_pair);
        

       pairsInterInter_ = [INT1(mask_to_pair) INT2(mask_to_pair)];
%        for j = 1:size(pairsInterInter_,1)
%             figure(21); 
%             hold off;
%             imshow(img);
%            hold on;
%            plot(intersections.coordinates2D(pairsInterInter_(j,:),1),...
%                 intersections.coordinates2D(pairsInterInter_(j,:),2), '*-');
%            plot([strokes_topology(intersections.strokes_indices(pairsInterInter_(j,1),1)).points2D(:).x],...
%                  [strokes_topology(intersections.strokes_indices(pairsInterInter_(j,1),1)).points2D(:).y]);
%            plot([strokes_topology(intersections.strokes_indices(pairsInterInter_(j,1),2)).points2D(:).x],...
%                  [strokes_topology(intersections.strokes_indices(pairsInterInter_(j,1),2)).points2D(:).y]);
%              
%            plot([strokes_topology(intersections.strokes_indices(pairsInterInter_(j,2),1)).points2D(:).x],...
%                  [strokes_topology(intersections.strokes_indices(pairsInterInter_(j,2),1)).points2D(:).y]);
%            plot([strokes_topology(intersections.strokes_indices(pairsInterInter_(j,2),2)).points2D(:).x],...
%                  [strokes_topology(intersections.strokes_indices(pairsInterInter_(j,2),2)).points2D(:).y]);
%        end
       
        pairsInterInter((end+1):(end+sum(mask_to_pair)),:) = [INT1(mask_to_pair) INT2(mask_to_pair)];        
    end

    
    
    
    
    
    
    
    
%% Method 2    
%     num_intersections = size(intersections.coordinates2D,1);
%     [INT1, INT2] = meshgrid(1:num_intersections);
%     
%     %Compute disatnces only once:
%     A = true(num_intersections);
%     A = triu(A,1);
%     INT1 = INT1(A);
%     INT2 = INT2(A);
% 
%     distances = sqrt(sum((intersections.coordinates2D(INT1,:) - intersections.coordinates2D(INT2,:)).^2,2));
% 
%     mask_to_pair =  (distances < intersections.accuracy_radius(INT1)) | ...
%                     (distances < intersections.accuracy_radius(INT2));
%     
%     distances = distances(mask_to_pair);
%     pairsInterInter = [INT1(mask_to_pair) INT2(mask_to_pair)];

%% Method 3
%     pairsInterInter = [];
%     distances = [];
%     for i1 = 1:num_intersections
%         inds_2 = (i1+1):num_intersections;
%         [INT1, INT2] = meshgrid(i1, inds_2);
%         distances_ = sqrt(...
%                         sum((intersections.coordinates2D(INT1,:) - ...
%                               intersections.coordinates2D(INT2,:)).^2, ...
%                           2)...
%                         );
%          mask_to_pair =  (distances_ < intersections.accuracy_radius(INT1)) | ...
%                     (distances_ < intersections.accuracy_radius(INT2));
%     
%         distances = [distances; distances_(mask_to_pair)];
%         pairsInterInter_ = [INT1(mask_to_pair) INT2(mask_to_pair)];
%         pairsInterInter = [pairsInterInter; pairsInterInter_];
%         
%     end

end
