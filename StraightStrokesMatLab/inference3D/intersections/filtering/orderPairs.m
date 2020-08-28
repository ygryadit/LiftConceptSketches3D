function pairs = orderPairs(pairs)
    ind_swap = find(pairs(:,2)<pairs(:,1));
    temp = pairs(ind_swap ,1);
    pairs(ind_swap ,1) = pairs(ind_swap ,2);
    pairs(ind_swap ,2) = temp;
    pairs = unique(pairs,'rows');
end