function direction_prior = getDirectionVec(line_grp)
    if line_grp == 1
       direction_prior = [1,0,0];
    elseif line_grp == 2
       direction_prior = [0,1,0];
    elseif line_grp == 3
       direction_prior = [0,0,1];     
    else
       error('Unknown direction, should be between 1 and 3');
    end
end