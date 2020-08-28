function reproject3Dto2D(img, cam_param, strokes_topology, intersections, fig_num, color)
close(figure(fig_num));
figure(fig_num); 
hold on;
% if exist('img', 'var')
% figure(fig_num); 
% hold off;
% imshow(img);
% hold on;
% plot(cam_param.principal_point(1),cam_param.principal_point(2),'*');
% % plot(center_2D(1),center_2D(2),'o');
% end

for i = 1:length(strokes_topology)
    if (strokes_topology(i).depth_assigned)
       points2D = zeros(size(strokes_topology(i).points3D,1),2);

        for j = 1:size(strokes_topology(i).points3D,1)
            try
                p2D = cam_param.P*[strokes_topology(i).points3D(j,:)'; 1.0];
            catch
                disp('') 
            end
            points2D(j,:) = p2D(1:2)./p2D(3);
        end
        plot(points2D(:,1),points2D(:,2),color);
    end
    
end


for i = 1:length(intersections)
    if isnan(intersections(i).coordinates3D)
        continue;
    end
    
    p2D = cam_param.P*[intersections(i).coordinates3D'; 1.0];

    p2D = p2D(1:2)./p2D(3);
            
    plot(p2D(1),p2D(2), '*');
end




end