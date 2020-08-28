function [ strokes_topology_out, intersections, cam_param] = readReconstructionJsonFelix( filepath )
   
   fid = fopen(filepath);
   raw = fread(fid,inf);
   str = char(raw');
   fclose(fid);
   object = jsondecode(str);
   
   
    strokes_topology = object.strokes_topology;
    intersections = object.intersections;
    cam_param = object.cam_param;

    
%     global ZUP;
%     if signz == -1
%         ZUP = false;
%     else
%         ZUP = true;
%     end

    signz = 1.0;
    
    
    for i = 1:length(strokes_topology)
        if ~isfield(strokes_topology, 'points3D')
            stroke = strokes_topology{i};
        else
            stroke = strokes_topology(i);
        end
        stroke.points3D = stroke.points3D*signz;
        stroke.primitive_geom = stroke.primitive_geom';                
        stroke.direction_vec = stroke.direction_vec';
        stroke.primitive_geom_3D = stroke.primitive_geom_3D;
        
        if isfield(stroke, 'candidate_lines') & ~isempty(stroke.candidate_lines)
            for j = 1:length(stroke.candidate_lines)
                stroke.candidate_lines(j).origin = stroke.candidate_lines(j).origin';
                stroke.candidate_lines(j).dir = stroke.candidate_lines(j).dir';
                stroke.candidate_lines(j).coordinates3D_prior = stroke.candidate_lines(j).coordinates3D_prior'*signz;
                for k = 1:length( stroke.candidate_lines(j).configurations)
                     config = stroke.candidate_lines(j).configurations(k);
                     config.inds_intrsctns = config.inds_intrsctns';
                     config.inds_intrsctns__assigned = config.inds_intrsctns__assigned';
                     config.inds_intrsctns__mult_cnddts = config.inds_intrsctns__mult_cnddts';
                     config.inds_intrsctns__mult_cnddts_ind = config.inds_intrsctns__mult_cnddts_ind';
                     config.inds_jnts_strks = config.inds_jnts_strks';
                     config.p_intrsctns_dists = config.p_intrsctns_dists';
                     config.inds_intrsctns__assigned = config.inds_intrsctns__assigned';
                     
                     stroke.candidate_lines(j).configurations(k) = config ;
                end
            end
            
        end
       
        strokes_topology_out(i).points3D = stroke.points3D;
        strokes_topology_out(i).depth_assigned = stroke.depth_assigned;
        strokes_topology_out(i).mean_pressure = stroke.mean_pressure;
        strokes_topology_out(i).line_group = stroke.line_group;
        strokes_topology_out(i).primitive_type = stroke.primitive_type;
    end

    for i = 1:length(intersections)
        intersections(i).coordinates3D = intersections(i).coordinates3D';
        intersections(i).coordinates2D = intersections(i).coordinates2D';
        intersections(i).coordinates3D = intersections(i).coordinates3D*signz;
        intersections(i).strokes_indices = intersections(i).strokes_indices';
        if ~isfield(strokes_topology(1), 'cnddts3D') 
            continue;
        end
        
        for j = 1:length(intersections(i).cnddts3D)
            if ~iscell(intersections(i).cnddts3D(j).cnddt_lns)
                intersections(i).cnddts3D(j).cnddt_lns = num2cell(intersections(i).cnddts3D(j).cnddt_lns);
            end
            if ~iscell(intersections(i).cnddts3D(j).cnfgrtns) 
                intersections(i).cnddts3D(j).cnfgrtns = num2cell(intersections(i).cnddts3D(j).cnfgrtns,3);
            end
            for k = 1:length(intersections(i).cnddts3D(j).cnfgrtns)
                if ~iscell(intersections(i).cnddts3D(j).cnfgrtns{k}) 
                 intersections(i).cnddts3D(j).cnfgrtns{k} = ...
                     num2cell(intersections(i).cnddts3D(j).cnfgrtns{k}(:)',2);
                end
            end
        end
        
    end

    
    
end


