% folder_files = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\_400';
% designer = 'student8';
% object = 'house';
% view = 'view1';
% intersections = assignAllIntersections(folder_files, ...
%                                                 designer,...
%                                                 object,...
%                                                 view)
global designer;
global object_name;

folder_files = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\_400';
designer = 'Professional5';
object_name = 'bumps';
view = 'view1';

[strokes_topology, intersections] = assignAllIntersections(folder_files, ...
                                                designer,...
                                                object_name,...
                                                view);


folder_files = 'C:\Users\yulia\Research\DesignSketch3D\Results_Cluster\v49\v49planes_400\_400';
designer = 'Professional6';
object_name = 'shampoo_bottle';
view = 'view1';
intersections = assignAllIntersections(folder_files, ...
                                                designer,...
                                                object_name,...
                                                view);
