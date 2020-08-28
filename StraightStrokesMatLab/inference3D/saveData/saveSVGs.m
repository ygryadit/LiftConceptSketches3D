function saveSVGs(strokes_topology)
    pressure_max = max(cat(1, strokes_topology(:).mean_pressure));
    thr_max_pressure = 1.0;
    if pressure_max < thr_max_pressure 
        scale_f = thr_max_pressure/pressure_max;
        for  i = 1:length(strokes_topology)
            strokes_topology(i).mean_pressure = strokes_topology(i).mean_pressure*scale_f;
        end
    end
    
    save_as_svg(strokes_topology)
    save_as_svg_straigt_strokes_vp_colorcoded(strokes_topology);
end

