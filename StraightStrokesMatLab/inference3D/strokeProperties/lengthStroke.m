function length = lengthStroke(stroke_points)
    length = 0.0;
    
    for p = 1:(size(stroke_points,1)-1)
        delta = stroke_points(p+1,:) - stroke_points(p,:);
        length = length + norm(delta);
    end

end