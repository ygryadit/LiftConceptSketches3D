function candidate_lines = updateListJointStrokes(candidate_lines)
    for i =1:length(candidate_lines)
        
        mask = arrayfun(@(x) ~isempty(x.list_jnt_strks), candidate_lines(i).configurations);
        
        try
        candidate_lines(i).list_jnt_strks = unique(cat(2, candidate_lines(i).configurations(mask).list_jnt_strks));
        catch
           disp(''); 
        end
    end
end



