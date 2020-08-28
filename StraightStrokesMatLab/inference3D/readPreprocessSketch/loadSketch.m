function [sketch,img,accuracy_radiuses] = loadSketch()

global filepath_sketch_json;
global filepath_sketch_img;
global sketch_height;


%% Read sketch json and image
sketch = readSketchJson(filepath_sketch_json);
sketch = keepOnlyNonRemovedStrokes(sketch);
sketch_height = sketch.canvas.width;
sketch.canvas.heigth = sketch.canvas.width; 
if isnan(sketch_height)
    sketch_height = findMaxXY(sketch);
    sketch.canvas.height = sketch_height;
    sketch.canvas.width = sketch_height;
  
end


img = readSketchImg(filepath_sketch_img, true);


%% Preprocess strokes:
%  Preprocess strokes to remove hooks and obtain simplified sampling.

DGP_THR = 0.5; %1.5/4.0;
for i = 1:length(sketch.strokes)
    sketch.strokes(i).mean_pressure = mean(cat(1,sketch.strokes(i).points(:).p));
    sketch.strokes(i).primitive_type = [];
end

[sketch.strokes] = computeAverageStrokesSpeed(sketch.strokes);

for i = 1:length(sketch.strokes)
    sketch.strokes(i).points = trimStraightStrokeCoordinatesDrawingField(sketch.strokes(i).points);
end

sketch.strokes = resampleAllStrokes(sketch.strokes, DGP_THR, img);
sketch.strokes = cleanSubdivideStrokes(sketch.strokes, img); 

%% Scale and center the sketch:

[sketch, img, scale] = scaleCenterSketch(sketch, img);
filepath_sketch_img = strrep(filepath_sketch_img, 'view1_concept_opaque.png', 'view1_centered.png');
imwrite(img, filepath_sketch_img);

%% Compute accuracy radiuses based on speed of the stroke:
%     [strokes_topology, mask_lines_kept]  = computeAverageStrokesSpeed(strokes_topology);
%     probabilitiesVP = probabilitiesVP(mask_lines_kept,:);
[accuracy_radiuses] = computeMergeThresholdStrokes(cat(1,sketch.strokes(:).speed), scale);


    
end

