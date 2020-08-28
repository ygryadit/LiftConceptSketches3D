function [ind_similar_cand_lines, inds_strs_same_intrsctns, candidate_lines] = ...
                findIndicesCandidateLinesSameIntersectionsAndDirections(...
                    candidate_lines,...
                    inds_strs_same_intrsctns,...
                    ind_cur_line,...
                    inds_strs_same_intrsctns_rest)
                   
    if isempty(inds_strs_same_intrsctns_rest)
        ind_similar_cand_lines = ind_cur_line;
        inds_strs_same_intrsctns = inds_strs_same_intrsctns_rest;
        return;
    end

    
    % Direction of the curent line:
    dir_cur_line = candidate_lines(ind_cur_line).dir;
    
    
    %Lines with similar directions
    dirs_rest_lines = cat(1,candidate_lines(inds_strs_same_intrsctns_rest).dir);

    dir_cur_line = repmat(dir_cur_line, size(dirs_rest_lines,1),1);

    cosines_lines = dot(dir_cur_line,dirs_rest_lines,2);

    % Unify the directions:
    ind_change_sign_dir = inds_strs_same_intrsctns_rest(cosines_lines < 0);

    for ll = ind_change_sign_dir'                
        candidate_lines(ll).dir = -candidate_lines(ll).dir;                
    end

    
%     ind_similar_cand_lines = inds_strs_same_intrsctns_rest(abs(cosines_lines) > 0.95);
    ind_similar_cand_lines = inds_strs_same_intrsctns_rest(abs(cosines_lines) > 0.99); %~ 8 degree
    ind_similar_cand_lines = [ind_cur_line; ind_similar_cand_lines];
    inds_strs_same_intrsctns = setdiff(inds_strs_same_intrsctns, ind_similar_cand_lines);
    
        
end
