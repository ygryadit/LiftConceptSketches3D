function [  p,...
            vp,...  % 3 dominant vanishing point
            vps_selected,... % additional vanishing points for sets of paralel lines
            inds_axis,...    % axis that is orthogonal to the line
            inds_active_lines ... % indices of lines that belond to a certain cluster]
            ] = getVanishingPoints(strokes_topology)
        
    global folderVP;
    global DEBUG;
    global fid;
    
    filenameVP = fullfile(folderVP, 'vps.mat');
    
    if DEBUG
        if exist(filenameVP, 'file')
            delete(filenameVP);
        end
    end
    
    if ~exist(filenameVP, 'file')
        fprintf(fid, '\t compute...\n');
        
        mask_lines = cat(1, strokes_topology(:).primitive_type) == 0;
        lines_coordinates2D = cat(1,strokes_topology(mask_lines).primitive_geom);
        
        
           [failed,...
            vp,...
            p,... % the probability to coverage towards one of three vanishing points or to neither
            vps_selected,... % additional vanishing points for sets of paralel lines
            inds_axis,...    % axis that is orthogonal to the line
            inds_active_lines ... % indices of lines that belond to a certain cluster
            ] = getVPSketch(lines_coordinates2D);           


        save(filenameVP, ...
                'failed', ...
                'vp',...
                'p',...
                'vps_selected',...
                'inds_axis',...
                'inds_active_lines');                                    
    else
        fprintf(fid, '\t load...\n');
        load(filenameVP,'vp',...
                'p', ...
                'vps_selected',...
                'inds_axis',...
                'inds_active_lines'); 
%         load(filenameVP,'failed','vp','p','VoteArrTemp'); 
%         folderpath = getFolderpathVP(designer, object_name, view);
%         visualizeVP(true,folderpath,img,vp,p,lines.coordinates2D);
    end
    vp = reshape(vp, 2, [])';
    
    
end