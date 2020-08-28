function coord = interpolate_x(seg, val_new)
    y2 = seg(2,2); 
    y1 = seg(1,2);
    
    x2 = seg(2,1);
    x1 = seg(1,1);
    
    dy = y2 - y1;
    y_grows = dy > 0;
    dx = x2 - x1;
    
    r = dy/dx;
    
    if y_grows
        coord = y1 + r*(val_new - x1);
    else
        coord = y2 + r*(val_new - x2);
    end

end