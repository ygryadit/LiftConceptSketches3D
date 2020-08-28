function [vp_selected,...
    ind_axis,...
    inds_inactive_lines,...
    inds_active_lines] = computeOneRotatedSetOfParallelLines(vp,lines_active,w,h,img)

    vp_selected = [];
    ind_axis = [];
    inds_inactive_lines = [];
    inds_active_lines = [];
    
    global SHOW_FIGS_PREPROCESS;
    %% Find all valid intersections:
    Xpnts = computeIntersectionPoints(lines_active);
    
    inds = ~isnan(Xpnts(:,1)) & ~isnan(Xpnts(:,2)) & ...
        ~isinf(Xpnts(:,1)) & ~isinf(Xpnts(:,2));
    Xpnts = Xpnts(inds,:);
    
    [Xpnts] = removeRedundantPoints(Xpnts,w,h);
    
    if SHOW_FIGS_PREPROCESS
        
        figure(4);
        imshow(img);
        hold on;
        plot(lines_active(:, [1 2])', lines_active(:, [3 4])','r');
        plot(Xpnts(:,1), Xpnts(:,2), 'ro');
        axis equal;
        
        plot(vp([1,2],1), vp([1,2],2),'b');
        plot(vp([3,2],1), vp([3,2],2),'b');
        plot(vp([1,3],1), vp([1,3],2),'b');
    end
    
    %% Keep only lines near lines connecting three main vanising points.
    
    distances = pointLineDistance3D(Xpnts, vp(1,:), vp(2,:));
    mask = distances < 0.2*w;
    
    distances = pointLineDistance3D(Xpnts, vp(2,:), vp(3,:));
    mask = mask | (distances < 0.2*w);
    
    distances = pointLineDistance3D(Xpnts, vp(1,:), vp(3,:));
    mask = mask | (distances < 0.2*w);
    Xpnts = Xpnts(mask,:);
    if isempty(Xpnts)
        return;
    end
    if SHOW_FIGS_PREPROCESS
        
        figure(5);
        imshow(img);
        hold on;
        plot(lines_active(:, [1 2])', lines_active(:, [3 4])','r');
        plot(Xpnts(:,1), Xpnts(:,2), 'ro');
        axis equal;
        
        plot(vp([1,2],1), vp([1,2],2),'b');
        plot(vp([3,2],1), vp([3,2],2),'b');
        plot(vp([1,3],1), vp([1,3],2),'b');
    end
    
    
    %% Compute votes:
%     VoteArr = computeLinesPointsVotes([lines_active;All_lines(~grp4,:)],Xpnts);
    VoteArr = computeLinesPointsVotes([lines_active],Xpnts);
    Vote=sum(VoteArr(1:size(lines_active,1),:),1);
      
    %Sort vanishing points and their votes:
    [vv, ii]=sort(Vote,'descend');
    Xpnts=Xpnts(ii,:);
    VoteArr = VoteArr(:,ii);
    
    vp_selected(1:2)=Xpnts(1,1:2);
    lines_votes_vp1 = VoteArr(:,1);
    
    max_length = max(lines_active(:,7));
    inds_active_lines    = find((lines_votes_vp1*max_length./lines_active(:,7))<0.8);
    inds_inactive_lines  = find((lines_votes_vp1*max_length./lines_active(:,7))>=0.8);    
    
    
     if SHOW_FIGS_PREPROCESS
       plot(vp_selected(1), vp_selected(2), '*g');
       plot(lines_active(inds_inactive_lines, [1 2])', lines_active(inds_inactive_lines, [3 4])','b');
     end
    
     ind_axis = findOrthPrior(vp_selected, vp);
end

function ind = findOrthPrior(vp_selected, vp)
  distance = zeros(3,1);
  distance(3) = pointLineDistance3D(vp_selected, vp(1,:), vp(2,:));
  
  distance(1) = pointLineDistance3D(vp_selected, vp(2,:), vp(3,:));    
  
  distance(2) = pointLineDistance3D(vp_selected, vp(1,:), vp(3,:));
    
  [~, ind] = min(distance);
end