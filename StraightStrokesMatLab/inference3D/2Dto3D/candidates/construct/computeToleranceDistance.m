function tolerance = computeToleranceDistance(strokes_topology, candidate_line)
%     based on approach from Lipson

%find distance to each stroke as an avarage disatance to it's endpoints,
%including half of the stroke itself, the minimum value is an accuracy
%radius.



lines_assigned = cat(1,strokes_topology([strokes_topology(:).depth_assigned]).primitive_geom_3D);
lines_assigned= lines_assigned';
lines_assigned = reshape(lines_assigned, 6, [])';

inds_mult_clns = find([strokes_topology(:).num_candidate_lines] > 0);

lines_non_assigned = [];

if isfield(strokes_topology, 'candidate_lines')
    clns = cat(2, strokes_topology(inds_mult_clns).candidate_lines);
    if ~isempty(clns)
        lines_non_assigned = cat(1, clns(:).coordinates3D_prior);
    end
end


line_cur = candidate_line.coordinates3D_prior;


lines = [lines_assigned; lines_non_assigned; line_cur];

distances = zeros(size(lines,1),2);
line_cur = reshape(line_cur', 3, 2)';
for i = 1:size(lines,1)
    
    distances(i,1)= point_seg_dist(line_cur(1,:),lines(i,1:3), lines(i,2:4));
    distances(i,2)= point_seg_dist(line_cur(2,:),lines(i,1:3), lines(i,2:4));
end
distances = mean(distances, 2);

tolerance = min(distances);
end

function d = point_seg_dist(P, A, B)
% Returns distance from a 3D point to a line segment
%
% Inputs
%
%   P : vector
%       Position in space
%
%   A : vector
%       Position of line segment endpoint
%
%   B : vector
%       Position of other line segment endpoint
%
% Outputs
%
%   d : double
%       Distance from point to line segment

AP = P - A; % Vector from A to point
AB = B - A; % Vector from A to B

% Project point onto line
P_line = A + dot(AP, AB) / dot(AB, AB) * AB;

if all(A < P_line) && all(P_line < B)
    % The point projected onto the line is in between A and B 

    % Projection of point onto segment is the same 
    % as projection of point onto line
    P_seg = P_line; 
else
    % The point projected onto the line is outside of A and B

    if all(P_line <= A)
        % The projected point is closer to A
        P_seg = A;  
    else
        % The projected point is closer to B
        P_seg = B;
    end
end

d = norm(P - P_seg); % Distance to line segment

end