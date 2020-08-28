function strokes_numbers = read_files(folder, filename)
   
    files = dir(fullfile(folder, filename));
    files = files(~[files(:).isdir]);
    files = {files.name};
    file = files{1};
    
    fid = fopen(fullfile(folder, file));
    tline = fgetl(fid);
    strokes_numbers = [];
    while ischar(tline)
%         disp(tline)
        C = cell2mat(textscan(tline,'o [%d]*'));
        if ~isempty(C)
            strokes_numbers(end+1) = C;
        end
        tline = fgetl(fid);
    end
    fclose(fid);    
end