function areContinious_ = areContinious(strokes_topology, ...
                                        ind_strk_1,...
                                        ind_strk_2)
        
    indcs_intrsctng_strks = strokes_topology(ind_strk_1).indcs_intrsctng_strks(strokes_topology(ind_strk_1).indcs_intrsctng_strks > ind_strk_1);
    indcs_intrsctng_strks = sort(indcs_intrsctng_strks, 'ascend');
    if ~isempty(indcs_intrsctng_strks)
        n = min(2, length(indcs_intrsctng_strks));
        areContinious_ = ismember(ind_strk_2, indcs_intrsctng_strks(1:n));
    else
        areContinious_  =false;
    end
                

end