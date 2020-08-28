function [ sketch ] = readSketchJson( filepath )
%    global fid;
%    fprintf(fid, 'readSketchJson(): read file: %s\n', filepath);   
   fidf = fopen(filepath);
   raw = fread(fidf,inf);
   str = char(raw');
   fclose(fidf);
   sketch = jsondecode(str);
%    fprintf(fid, 'done \n');
end


