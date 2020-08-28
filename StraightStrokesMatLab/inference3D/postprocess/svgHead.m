function fidSVG = svgHead(fidSVG)
global sketch_height;
str = sprintf('<?xml version="1.0" encoding="utf-8" ?>\n'); 
fwrite(fidSVG, str);
str = sprintf('<svg baseProfile="full" height="%d" version="1.1" viewBox="0,0,%d,%d" width="%d" xmlns="http://www.w3.org/2000/svg" xmlns:ev="http://www.w3.org/2001/xml-events" xmlns:xlink="http://www.w3.org/1999/xlink"><defs><style type="text/css"><![CDATA[\n', ...
    sketch_height,sketch_height,sketch_height,sketch_height);
fwrite(fidSVG, str);
str = sprintf('\t.background { fill: white; }\n');
fwrite(fidSVG, str);
str = sprintf('\t.line { stroke: firebrick; stroke-width: .1mm; }\n');
fwrite(fidSVG, str);
str = sprintf('\t.blacksquare { fill: indigo; }\n');
fwrite(fidSVG, str);
str = sprintf('\t.whitesquare { fill: white; }\n');
fwrite(fidSVG, str);
str = sprintf(']]></style></defs>\n');
fwrite(fidSVG, str);
end