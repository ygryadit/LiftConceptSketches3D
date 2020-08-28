function [x_vals,  y_vals] = ...
            updateOneClipCoordinateMax(x_vals,...
                                    y_vals,...
                                    inds_out,...
                                    grows_dim,...
                                    interp_coord,...
                                    max_b)
                                
num_points = length(y_vals);
update = false;
if ~isempty(inds_out)
    if grows_dim
        vert_num = inds_out(1);

        if vert_num ~= 1
            update = true;
            seg = [x_vals(vert_num) y_vals(vert_num); ...
                   x_vals(vert_num-1) y_vals(vert_num-1)];  
        end
    else
        vert_num = inds_out(1);

        if vert_num ~= num_points
            update = true;
            seg = [x_vals(vert_num)  y_vals(vert_num); ...
                   x_vals(vert_num+1) y_vals(vert_num+1)];

        end
    end

    if  update
        if strcmp(interp_coord, 'y')
            x_vals(vert_num) = interpolate_y(seg, max_b);
            y_vals(vert_num) = max_b;
        else
            x_vals(vert_num) = max_b;
            y_vals(vert_num) = interpolate_x(seg, max_b);
        end
    end

end

end