function ind = findCoincidingLines(intersections_strokes_indices, lines, img)

global SHOW_FIGS_PREPROCESS

dir1 =  lines(intersections_strokes_indices(:,1), [2,4]) - ...
        lines(intersections_strokes_indices(:,1), [1,3]);
dir2 =  lines(intersections_strokes_indices(:,2), [2,4]) - ...
        lines(intersections_strokes_indices(:,2), [1,3]);    

dir1 = dir1./repmat(sqrt(sum(dir1.^2,2)), 1, 2);
dir2 = dir2./repmat(sqrt(sum(dir2.^2,2)), 1, 2);

cos_dirs = dot(dir1, dir2, 2);


% disp(acos(cos_dirs)/pi*180)
ind = find(abs(cos_dirs) > cos(pi/36));





if SHOW_FIGS_PREPROCESS
%     figure;
%     for i = ind'
% 
%         hold off;
%         imshow(img);
%         hold on;
%         si = intersections_strokes_indices(i,1);
%         sj = intersections_strokes_indices(i,2);
%         plot(lines(si,[1,2]), lines(si,[3,4]), ':', 'LineWidth', 4);
%         plot(lines(sj,[1,2]), lines(sj,[3,4]), ':', 'LineWidth', 4);
%     end
end

end