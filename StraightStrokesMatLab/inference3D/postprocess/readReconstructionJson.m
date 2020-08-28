function [ strokes_topology, intersections, cam_param] = readReconstructionJson( filepath )
   
   fid = fopen(filepath);
   raw = fread(fid,inf);
   str = char(raw');
   fclose(fid);
   object = jsondecode(str);
   
   
    strokes_topology = object.strokes_topology;
    intersections = object.intersections;
    cam_param = object.cam_param;
    signz = object.signz;
    
%     global ZUP;
%     if signz == -1
%         ZUP = false;
%     else
%         ZUP = true;
%     end

    signz = 1.0;
    
    for i = 1:length(strokes_topology)
        strokes_topology(i).points3D = strokes_topology(i).points3D*signz;
        if isfield(strokes_topology(i), 'primitive_geom')
            strokes_topology(i).primitive_geom = strokes_topology(i).primitive_geom';                
        end
        if isfield(strokes_topology(i), 'direction_vec')
            strokes_topology(i).direction_vec = strokes_topology(i).direction_vec';
        end
        if isfield(strokes_topology(i), 'primitive_geom_3D')
            strokes_topology(i).primitive_geom_3D = strokes_topology(i).primitive_geom_3D;
        end
        
        if isfield(strokes_topology(i), 'candidate_lines') & ~isempty(strokes_topology(i).candidate_lines)
            for j = 1:length(strokes_topology(i).candidate_lines)
                strokes_topology(i).candidate_lines(j).origin = strokes_topology(i).candidate_lines(j).origin';
                strokes_topology(i).candidate_lines(j).dir = strokes_topology(i).candidate_lines(j).dir';
                strokes_topology(i).candidate_lines(j).coordinates3D_prior = strokes_topology(i).candidate_lines(j).coordinates3D_prior'*signz;
                for k = 1:length( strokes_topology(i).candidate_lines(j).configurations)
                     config = strokes_topology(i).candidate_lines(j).configurations(k);
                     config.inds_intrsctns = config.inds_intrsctns';
                     config.inds_intrsctns__assigned = config.inds_intrsctns__assigned';
                     config.inds_intrsctns__mult_cnddts = config.inds_intrsctns__mult_cnddts';
                     config.inds_intrsctns__mult_cnddts_ind = config.inds_intrsctns__mult_cnddts_ind';
                     config.inds_jnts_strks = config.inds_jnts_strks';
                     config.p_intrsctns_dists = config.p_intrsctns_dists';
                     config.inds_intrsctns__assigned = config.inds_intrsctns__assigned';
                     
                     strokes_topology(i).candidate_lines(j).configurations(k) = config ;
                end
            end
            
        end
        
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


