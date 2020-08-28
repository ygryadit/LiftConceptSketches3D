function data = loadLabeling(filepath)
   
   fid = fopen(filepath);
   raw = fread(fid,inf);
   str = char(raw');
   fclose(fid);
   data = jsondecode(str);
   
end