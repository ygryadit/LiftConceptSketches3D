function visualizeVP(do_visualize,folderpath_vp, img,vp,p,All_lines)
    if ~do_visualize
        return;
    end
%     figure(2);
%     subplot(1,2,1);
%     plot(vp(1),vp(2),'r*');
%     hold on;
%     imagesc(img);hold on;
%     plot(vp(1),vp(2),'r*');
%     plot(vp(3),vp(4),'g*');
%     plot(vp(5),vp(6),'b*');

    [vv, linemem] = max(p,[],2);
        
    grp1=find(linemem==1);
    grp2=find(linemem==2);
    grp3=find(linemem==3);
    grp4=find(linemem==4);
%     plot(All_lines(grp1, [1 2])', All_lines(grp1, [3 4])','r');
%     plot(All_lines(grp2, [1 2])', All_lines(grp2, [3 4])','g');
%     plot(All_lines(grp3, [1 2])', All_lines(grp3, [3 4])','b');
%     plot(All_lines(grp4, [1 2])', All_lines(grp4, [3 4])','c');
%     axis ij;axis equal;
%     
%     subplot(1,2,2);
%     imagesc(img);hold on;
%     plot(All_lines(grp1, [1 2])', All_lines(grp1, [3 4])','r');
%     plot(All_lines(grp2, [1 2])', All_lines(grp2, [3 4])','g');
%     plot(All_lines(grp3, [1 2])', All_lines(grp3, [3 4])','b');
%     plot(All_lines(grp4, [1 2])', All_lines(grp4, [3 4])','c');
%     
%     axis ij;axis equal;
%     saveas(2, [folderpath_vp '.fig']);
    
%     close all;
    figure(10);
    plot([1 size(img,1) 1 size(img,1)], [1 1 size(img,2) size(img,2)], 'w');    
    imagesc(img);hold on;
    axis off; axis equal;
    plot(All_lines(grp4, [1 2])', All_lines(grp4, [3 4])','Color',  [255, 173, 51]/255.0, 'LineWidth', 0.5, 'LineStyle', ':');
    plot(All_lines(grp1, [1 2])', All_lines(grp1, [3 4])', 'Color', 'm', 'LineWidth', 0.5);
    plot(All_lines(grp2, [1 2])', All_lines(grp2, [3 4])', 'Color', [0, 255, 255]/255.0, 'LineWidth', 0.5);
    plot(All_lines(grp3, [1 2])', All_lines(grp3, [3 4])', 'Color', [27, 192, 27]/255.0, 'LineWidth', 0.5);
    
    
    
    clear 'centers'
    centers(1,:) = (All_lines(grp1, 1)' + All_lines(grp1, 2)')/2.0;
    centers(2,:) = (All_lines(grp1, 3)' + All_lines(grp1, 4)')/2.0;   
    [x,y] = scaleCoordiantes(centers, vp(1:2), size(img,1));
    plot([centers(1,:); x], [centers(2,:); y], 'Color', 'm', 'LineStyle', ':', 'LineWidth', 1.0);
    clear 'centers'
    
    centers(1,:) = (All_lines(grp2, 1)' + All_lines(grp2, 2)')/2.0;
    centers(2,:) = (All_lines(grp2, 3)' + All_lines(grp2, 4)')/2.0;   
    [x,y] = scaleCoordiantes(centers, vp(3:4), size(img,1));
    plot([centers(1,:); x], [centers(2,:); y], 'Color', [0, 255, 255]/255.0, 'LineStyle', ':', 'LineWidth', 1.0);
    
    clear 'centers'
    centers(1,:) = (All_lines(grp3, 1)' + All_lines(grp3, 2)')/2.0;
    centers(2,:) = (All_lines(grp3, 3)' + All_lines(grp3, 4)')/2.0;   
    [x,y] = scaleCoordiantes(centers, vp(5:6), size(img,1));
    plot([centers(1,:); x], [centers(2,:); y], 'Color', [27, 192, 27]/255.0, 'LineStyle', ':', 'LineWidth', 1.0);
    
    
    saveas(10, fullfile(folderpath_vp, 'lines_convergence.png'));
    hold off;
end


function [x,y] = scaleCoordiantes(centers, vp, width)
    
    if (vp(1) > width)
        x = width*ones(1, length(centers(1,:)));              
        y = vp(2) + ( centers(2,:)-vp(2)).*(vp(1)-width)./(vp(1) - centers(1,:));
    elseif (vp(1) < 0)
        x = zeros(1, length(centers(1,:)));
        y = vp(2) + (centers(2,:) - vp(2)).*abs(vp(1))./(centers(1,:) + abs(vp(1)));
    elseif (vp(2) > width)
        y = width*ones(1, length(centers(1,:)));              
        x = vp(1) + ( centers(1,:)-vp(1)).*(vp(2)-width)./(vp(2) - centers(2,:));
    elseif (vp(2) < 0)
        y = zeros(1, length(centers(1,:)));
        x = vp(1) + (centers(1,:) - vp(1)).*abs(vp(2))./(centers(2,:) + abs(vp(2)));
        
    else
        
        
        x = vp(1)*ones(1, length(centers(1,:)));
        y = vp(2)*ones(1, length(centers(1,:)));
    end    
end    
