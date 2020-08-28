function lines_pairs = rearrangePairs(lines_pairs)
   lines_pairs = unique(lines_pairs, 'rows');
   mask_diff = find(lines_pairs(:,1)~=lines_pairs(:,2));
   lines_pairs = lines_pairs(mask_diff,:);
   
   inds = find(lines_pairs(:,2) < lines_pairs(:,1));
   temp = lines_pairs(inds,1);
   lines_pairs(inds,1) = lines_pairs(inds,2);
   lines_pairs(inds,2) = temp;
    
end