function saveMP4s(strokes_topology)
    global folder_save;
    global SAVE_PROCESS;
    global SHOW_FIGS;
    global SAVE_SVGS;
    global frame_count;
    
    if SAVE_PROCESS && SHOW_FIGS
        frame_count = frame_count +1;
        f1 =  fullfile(folder_save, '3D_decision');
        f2 =  fullfile(folder_save, '3D_candidates');
        f3 =  fullfile(folder_save, '2D');
        
        if ~exist(f1, 'dir')
            mkdir(f1);
        end
        if ~exist(f2, 'dir')
            mkdir(f2);
        end
        if ~exist(f3, 'dir')
            mkdir(f3);
        end
        
        filename = sprintf('%05d.png',frame_count);
        
        fl_3D_decision        = fullfile(f1,filename);
        fl_3D_candidates      = fullfile(f2,filename);
        fl_2D                 = fullfile(f3,filename);
        %% ============= Save to gif 3D decision ===================================
%         filetemp = fullfile(folder_save,'temp.mp4');
        writeFrames(9, fl_3D_decision);   
        writeFrames(11, fl_3D_candidates);
        writeFrames(2, fl_2D); 
    end
    
    if SAVE_SVGS
        saveSVGs(strokes_topology);
    end
end


