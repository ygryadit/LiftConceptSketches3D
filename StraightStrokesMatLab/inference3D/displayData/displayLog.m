function displayLog(message)
    
    global DISPLAY_INFO;
    global fid;
    if DISPLAY_INFO
        fprintf(fid, [message '\n']);
    end    
end