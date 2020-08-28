% function point_3D = estimateOneCoordinate(ind_coordinate, node_coordinates, node_prev_3D, P)
% 
% Input:
%   ind_coordinate:
%       possible values are 1,2 and 3. Defines which coordiante of the 
%       (x,y,z) vector should be estiamted.   
%   node_coordinates:
%       2D corrdiantes of the node        
%   node_prev_3D:
%       the 3D coordiantes of the node lying on the line with
%       'ind_coordinate' direction.
%   P:
%       projection matrix.

function point_3D = estimateOneCoordinate(ind_coordinate, node_coordinates, node_prev_3D, P, img)
    point_3D = node_prev_3D;
    u = node_coordinates(1);
    v = node_coordinates(2);
    ind_other = setxor(1:3, ind_coordinate);
    
%     point_3D(ind_coordinate) = 1.0/( P(1, ind_coordinate) - P(3, ind_coordinate)*u) * ...
%             ( u*(dot(P(3,ind_other),node_prev_3D(ind_other)) + P(3,4)) - ...
%               P(1,4) - dot(P(1,ind_other),node_prev_3D(ind_other))...
%             );
%     pp = P*[point_3D 1]';
%     pp = pp./pp(3)
% 
%     point_3D(ind_coordinate) = 1.0/( P(2, ind_coordinate) - P(3, ind_coordinate)*v) * ...
%             ( v*(dot(P(3,ind_other),node_prev_3D(ind_other)) + P(3,4)) - ...
%               P(2,4) - dot(P(2,ind_other), node_prev_3D(ind_other))...
%             );
%         
%     pp = P*[point_3D 1]';
%     pp = pp./pp(3)
    
    A = [(  P(1, ind_coordinate) - P(3, ind_coordinate)*u);...
            P(2, ind_coordinate) - P(3, ind_coordinate)*v];
    b = [u*(dot(P(3,ind_other),node_prev_3D(ind_other)) + P(3,4)) - ...
               P(1,4) - dot(P(1,ind_other),node_prev_3D(ind_other));
         v*(dot(P(3,ind_other),node_prev_3D(ind_other)) + P(3,4)) - ...
               P(2,4) - dot(P(2,ind_other), node_prev_3D(ind_other))];

    point_3D(ind_coordinate) = A\b;
    
%%  Visualisation:    
%     if exist('img', 'var')
%      color = ['r', 'g', 'b'];
%      pp = P*[point_3D 1]';
%      pp = pp./pp(3);
%      figure(2); 
%      imshow(img); hold on;
%      plot(pp(1), pp(2), [color(ind_coordinate) 'o']);
%      hold on;
%      plot(u, v, [color(ind_coordinate) '*']);
%      hold off;
%     end
end

