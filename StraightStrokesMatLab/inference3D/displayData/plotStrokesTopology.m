function plotStrokesTopology(fig_num, strokes_assigned) 
%     figure(fig_num);
%     holf off;
    num_strokes = length(strokes_assigned);
    
    colors = uint8(colormap(parula(num_strokes))*255);
   
    h = figure(fig_num);
     
%     set(h,'Position', [1168 74 1386 1195]);
%     set(h,'Position', [3522         178         958        1074]);
%     hold off;
   
   
%         view(-30.2000, 21.2000);
    set(h, 'Name', '3D candidate intersections');
   

    
    for i = 1:length(strokes_assigned)
        plot3(  strokes_assigned(i).points3D(:,1),...
                strokes_assigned(i).points3D(:,2),...
                strokes_assigned(i).points3D(:,3), ...
                'LineWidth', 1.0,...
                'Color', colors(i, :));
         hold on;
    end
    grid on;
    hold on;
    h = colorbar;
     
    b = num2str([1:num_strokes]');
    c = cellstr(b);
    h.Ticks = [0:(num_strokes-1)]/(num_strokes-1);
    h.TickLabels = c;
    
     axis equal;
    global ZUP
    if ~ZUP
        axis ij
        set(gca, 'Zdir', 'reverse')
    end
    global CameraView;
    view(3);
    set(gca, 'CameraPosition', [-1.6560 -3.3123 -1.6620]);
    set(gca, 'CameraTarget', [-0.1179 0.0548 0.1537]);
    
%     view(CameraView(1),CameraView(2));
%     camva(11.2020);
%     campos([-1.0327   -0.7573   -0.5696]);
%     camtarget([0.0319   0.0426    0.1172]);
%     camup([0     0    -1]);
end