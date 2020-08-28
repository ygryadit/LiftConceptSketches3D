% function P_plane = checkLinePlanarity(candidate_line)
% 
% Description: 
% each configuration stores all the normals of the planes formed by two
% strokes intersecting at each intersection that belong to the stroke
% defined by a configuration. We then check if the current candidate line
% and its configuration lie in some plane that exists for the stroke it
% intersects.
% 
% Input:
%   candidate_line:
%       one candidate lines with configurations
% 
% Output:
%   P_plane:
%       a vector of probabilities for each of the configurations that the
%       current candidate line lies in the existing plane in one of the
%       intersections.

function P_plane = checkLinePlanarity(candidate_line, ...
                                      intersections, ...
                                      strokes_topology, ...
                                      cur_stroke_ind)

    num_confgrtns = length(candidate_line.configurations);
    P_plane = zeros(num_confgrtns, 1);

    for i = 1:num_confgrtns
        configuration = candidate_line.configurations(i);

        P_plane(i) =  checkPlanarityGivenConfiguration(candidate_line.dir, ... the stoke which is being updated
                                            configuration, ... 
                                            intersections,...
                                            cur_stroke_ind,... the stoke which is being updated
                                            strokes_topology);

    end

end

