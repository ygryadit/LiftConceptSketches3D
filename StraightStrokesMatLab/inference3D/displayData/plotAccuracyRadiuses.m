function plotAccuracyRadiuses(img, strokes)
    
    global folder_save_imgs;

    hf = figure; 
    imshow(img);
    hold on;
    colors = colorcube(length(strokes));    
    for i = 1:length(strokes)
        plot(cat(1,strokes(i).points2D.x),cat(1,strokes(i).points2D.y), 'Color', colors(i,:));
        plotCircle(strokes(i).points2D(1).x, strokes(i).points2D(1).y, strokes(i).accuracy_radius, colors(i,:));
        num_points = length(strokes(i).points2D);
        plotCircle(strokes(i).points2D(num_points).x, strokes(i).points2D(num_points).y, strokes(i).accuracy_radius, colors(i,:));
    end

    set(hf, 'Name', 'Accuracy radius');

    save_png_file = fullfile(folder_save_imgs, 'accuracy_radius.png');
    saveas(hf, save_png_file);
    
end