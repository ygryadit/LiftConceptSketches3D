function fov = focalLength2FieldOfView(f, width)
    fov = 2*atand(width/(2.0*f));
end

