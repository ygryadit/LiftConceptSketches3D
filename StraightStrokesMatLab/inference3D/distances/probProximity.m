function E_distance = probProximity(distances, accuracy)   
    E_distance = exp(-distances.^2./(2*accuracy.^2));
end
