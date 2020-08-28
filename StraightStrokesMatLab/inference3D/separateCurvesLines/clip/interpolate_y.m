function coord = interpolate_y(seg, val_new)
    y2 = seg(2,2); 
    y1 = seg(1,2);
    
    x2 = seg(2,1);
    x1 = seg(1,1);
    
    dx = x2 - x1;
    dy = y2 - y1;
    x_grows = dx > 0;
    
    
    r = dx/dy;
    
    if x_grows
        coord = x1 + r*(val_new - y1);
    else
        coord = x2 + r*(val_new - y2);
    end

end