function normals = selectNormalsExcludeCurStroke(normals, cur_strk_ind)
    intrsctng_stks = normals(:,4);
    intrsctng_stks_mask = find(intrsctng_stks ~= cur_strk_ind);        
    normals = normals(intrsctng_stks_mask, 1:3);
end