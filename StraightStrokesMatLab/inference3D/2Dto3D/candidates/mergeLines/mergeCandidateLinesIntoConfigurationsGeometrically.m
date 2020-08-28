function candidate_lines_merged = mergeCandidateLinesIntoConfigurationsGeometrically(candidate_lines_new, ...
                                      cur_stroke,...
                                      cam_param,...
                                      strokes_topology, ...
                                       intersections )

                                  
num_lines = length(candidate_lines_new);
if num_lines < 2
    candidate_lines_merged = candidate_lines_new;
    return;
end

[IND1, IND2] = meshgrid(1:num_lines);

I = ones(num_lines);
U = triu(I);
IND1 = IND1(:);
IND2 = IND2(:);
U = logical(U(:));
IND1 = IND1(U);
IND2 = IND2(U);

lns1 = cat(1,candidate_lines_new(IND1).length3D);
lns2 = cat(1,candidate_lines_new(IND2).length3D);

%% Threshold by lines length:
vals = min(lns1./lns2, lns2./lns1);
INDSPairs = find(vals > 0.85);

IND1 = IND1(INDSPairs);
IND2 = IND2(INDSPairs);

lns1 = cat(1,candidate_lines_new(IND1).length3D);
lns2 = cat(1,candidate_lines_new(IND2).length3D);

%% Threshold by lines distances:
thr = 0.01;
thrs1 = lns1*thr;
thrs2 = lns2*thr;

thrs = min([thrs1, thrs2],[],2);
diffrncs = abs(cat(1,candidate_lines_new(IND1).coordinates3D_prior) - ...
               cat(1,candidate_lines_new(IND2).coordinates3D_prior));

distncs = diffrncs.^2;
distncs = 0.5*(sqrt(sum(distncs(:,1:3),2)) + sqrt(sum(distncs(:,4:6),2)));

angular_sym = dot(cat(1,candidate_lines_new(IND1).dir), cat(1,candidate_lines_new(IND2).dir), 2);

inds_merge = find((distncs < thrs) & abs(angular_sym) > 0.99);
% inds_merge = find(abs(angular_sym) > 0.99);

inds_pairs = [IND1(inds_merge) IND2(inds_merge)];

% groups = {};
candidate_lines_merged = [];

while ~isempty(inds_pairs)
    inds_lns = inds_pairs(1,1);
    [inds_lns,inds_pairs] = find_paired_lines(inds_lns, inds_pairs);
%     groups{end+1} = inds_lns;
    
    [candidate_line_] = mergeIntoSingleLineWithMultipleConfirations(candidate_lines_new,...
                                        inds_lns,...
                                        cur_stroke,...
                                        cam_param);
    if isempty(candidate_lines_merged)
        candidate_lines_merged = candidate_line_;
    else
        candidate_lines_merged(end+1) = candidate_line_;
    end
end




end

function [inds_lns,inds_pairs] = find_paired_lines(inds_lns, inds_pairs)
    mask1 = ismember(inds_pairs(:,1), inds_lns);
    mask2 = ismember(inds_pairs(:,2), inds_lns);
    mask = mask1|mask2;
    if (sum(mask) == 0)
        return;
    end
    
    inds_lns = unique([inds_lns; inds_pairs(mask1,2)]);
    inds_lns = unique([inds_lns; inds_pairs(mask2,1)]);
    
    inds_pairs = inds_pairs(find(~mask), :);
    [inds_lns,inds_pairs] = find_paired_lines(inds_lns, inds_pairs);
end