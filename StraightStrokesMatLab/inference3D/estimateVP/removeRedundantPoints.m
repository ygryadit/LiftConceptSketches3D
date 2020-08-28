%Copyright (c) October,15 2008 by Varsha Hedau, UIUC.  All rights reserved.
% function [Xpnts2,Vote2,VoteArr2] = removeRedundantPoints(Xpnts,Vote,VoteArr,w,h)
function [Xpnts2] = removeRedundantPoints(Xpnts,w,h)
%% Remove points that are too close to image center
% fprintf('number of points initial: %d\n', size(Xpnts, 1));
inds = find(Xpnts(:,1) > 0.8*w | Xpnts(:,1) < 0.2*w |...
       Xpnts(:,2) > 0.8*h | Xpnts(:,2) < 0.2*h);

% fprintf('number of points initial: %d\n', size(inds, 1));

Xpnts   = Xpnts(inds, :);
% Vote    = Vote(inds, :);
% VoteArr = VoteArr(inds, :);

%% Remove redundant points
currid=1;

Xpnts2  = Xpnts;
% Vote2   = Vote;
% VoteArr2= VoteArr;



while size(Xpnts,1)>0
    Xpnts2(currid,:)    = Xpnts(1,:);
%     Vote2(currid)       = Vote(1);
%     VoteArr2(:,currid)  = VoteArr(:,1);
    
  
    
    dists = (Xpnts(1,1)-Xpnts(:,1)).^2 + ...
            (Xpnts(1,2)-Xpnts(:,2)).^2;
        
    if sqrt((Xpnts(1,1)-w)^2+(Xpnts(1,2)-h)^2)/sqrt((w/2)^2+(h/2)^2) < 1
        % Higher precision close to the sketch:
        thres=20;
    else
        % Precision thta varies as we go further from the center:
        thres=30*(sqrt((Xpnts(1,1)-w)^2+(Xpnts(1,2)-h)^2)/sqrt((w/2)^2+(h/2)^2));
    end
    

    inds = find(dists<thres);
    num_points = size(inds,1);
    Xpnts2(currid,1) = sum(Xpnts(inds,1))/num_points;
    Xpnts2(currid,2) = sum(Xpnts(inds,2))/num_points;

    inds = find(dists>thres);
    Xpnts   = Xpnts(inds,:);
    
    
%     Vote    = Vote(inds,:);
%     VoteArr = VoteArr(:,inds);
    
   currid=currid+1;
end

Xpnts2 = Xpnts2(1:currid-1,:);
% Vote2 = Vote2(1:currid-1);
% VoteArr2 = VoteArr2(:,1:currid-1);


end
