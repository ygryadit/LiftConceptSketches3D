function [cost, inds_intrsctns, p_intrsctns_distances] = ...
                computeCost2DCoverage(...
                                inds_intrsctns,...
                                p_intrsctns_distances, ...
                                intersections,...
                                primitive_geom,...
                                line_length2D)
    
    
    
%     inds_intrsctns  = inds_intrsctns(~cat(1,intersections(inds_intrsctns).collinear));    
    num_intersections = length(inds_intrsctns);
    
    if num_intersections <= 1
        cost = 0;
        return;
    end
    
    
    
    inter_coord2D = cat(1,intersections(inds_intrsctns).coordinates2D);

    % Sort 2D coordinates along the geometric prior:

    [inter_coord2D, ...
      ind_sorted] = getSorted2DIntersectionCoordiantesAlongLineDir(inter_coord2D, primitive_geom([1,3]), primitive_geom([2,4]));

    inds_intrsctns = inds_intrsctns(ind_sorted);
    p_intrsctns_distances = p_intrsctns_distances(ind_sorted);

    
%     inds_intrsctns_nc  = inds_intrsctns(~cat(1,intersections(inds_intrsctns).collinear));
%     p_intrsctns_distances_nc = p_intrsctns_distances(~cat(1,intersections(inds_intrsctns).collinear));
    inds_intrsctns_nc = inds_intrsctns;
    p_intrsctns_distances_nc = p_intrsctns_distances;

    num_intersections = length(inds_intrsctns_nc);
    
    if num_intersections <= 1
        cost = 0;
        return;
    end
    
    % Compute maxx coverage:
    pairs = nchoosek(1:num_intersections, 2);
    pairs_inter_indices = inds_intrsctns_nc(pairs);

        
    ratios = sqrt(sum(( cat(1,intersections(pairs_inter_indices(:,1)).coordinates2D) - ...
                        cat(1,intersections(pairs_inter_indices(:,2)).coordinates2D) ).^2,2))/line_length2D;


    pSegments = ratios.*...
               p_intrsctns_distances_nc(pairs(:,1))'.*...
               p_intrsctns_distances_nc(pairs(:,2))';



    cost = computeWeightedAverage(inds_intrsctns_nc, pairs_inter_indices, pSegments);

    % Remap:
    cost = mapCost(cost);        
end


function cost = mapCost(cost)
    if cost > 0.85
        cost = 1.0;
    else 
        cost = exp(-(0.85 - cost).^2./(2*0.25^2));
    end
end

function max_cost = computeWeightedAverage(inds_intrsctns, pairs_inter_indices, pSegments)
    ind_pair = find( (pairs_inter_indices(:,1) == inds_intrsctns(1) & pairs_inter_indices(:,2) ==  inds_intrsctns(end) ) | ...
                     (pairs_inter_indices(:,2) == inds_intrsctns(1) & pairs_inter_indices(:,1) ==  inds_intrsctns(end) ) );
                 
    
    max_cost = pSegments(ind_pair);
    if length(inds_intrsctns) == 2
%         disp(max_cost);
        return;
    else
        ind_pair = find( (pairs_inter_indices(:,1) == inds_intrsctns(end-1) & pairs_inter_indices(:,2) ==  inds_intrsctns(end) ) | ...
                     (pairs_inter_indices(:,2) == inds_intrsctns(end-1) & pairs_inter_indices(:,1) ==  inds_intrsctns(end) ) );
        
        cost_left = computeWeightedAverage(inds_intrsctns(1:end-1), pairs_inter_indices, pSegments);
        cost_right = pSegments(ind_pair);
    
        max_cost = max(max_cost, cost_left+cost_right);
        
%         disp(inds_intrsctns(1));
%         disp(inds_intrsctns(end));
%         disp(max_cost);
    end
end