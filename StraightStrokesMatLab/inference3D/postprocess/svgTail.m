function svgTail(fidSVG)
str = sprintf('</svg>');
fwrite(fidSVG, str);

fclose(fidSVG);

end