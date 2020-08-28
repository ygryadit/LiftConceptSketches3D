% function [vp, lines_votes_vp1, active_lines, inactive_lines] = getDominantVP(lines)
% 
% Slects the vanishing point with the majority of vting lines, Then it
% devides the lines to those that are likely to converge towards this line
% and not.
% 
% Output:
%       vp:     
%           1x2 2D coordinates of a vanishing point
%       lines_votes_vp1: 
%           votes of each of the line stroke for the selcted vanishing
%           point
%       active_lines:
%           strokes numbers that are not likely to converge towards current
%           vanishing point
%      inactive_lines: 
%           strokes numbers that are likely to converge towards current
%           vanishing point
% 
% Input:
%       lines:
%           all input lines

function [vp, lines_votes_vp1, active_lines, inactive_lines] = getDominantVP(lines)
    global fid;
    global DISPLAY_INFO;
    global SHOW_FIGS_PREPROCESS;
      global filepath_sketch_img;
    % Compute possible candidate lines:
    %     Xpnts:      an array of candidate vanishing 
    Xpnts = computeIntersectionPoints(lines);
    inds = find(~isnan(Xpnts(:,1)) & ~isnan(Xpnts(:,2)) & ...
                ~isinf(Xpnts(:,1)) & ~isinf(Xpnts(:,2)));
    Xpnts = Xpnts(inds,:);
          
    % Computing votes for every point from all All_lines:
    %       VoteArr:    num_lines by num_vp matrix
    VoteArr = computeLinesPointsVotes(lines,Xpnts);
    Vote = sum(VoteArr,1); % a summary vote for each vanishing point
    
    % Get the first point & remove the All_lines of this point:
    %       ii:         sorted votes for each vanishing point
    [~, ii]=sort(Vote,'descend');    
    vp_selected(1:2)=Xpnts(ii(1),1:2);
    lines_votes_vp1 = VoteArr(:,ii(1));
    
    % Todo: remove scaling for max length
    max_length = max(lines(:,7));
    active_lines    = find((lines_votes_vp1*max_length./lines(:,7))<0.8);
    inactive_lines  = find((lines_votes_vp1*max_length./lines(:,7))>=0.8);    
    
%     if SHOW_FIGS_PREPROCESS
%         figure;
%         img = readSketchImg(filepath_sketch_img);
%         imshow(img);
%         hold on;        
%         plot(lines(inactive_lines, [1 2])', lines(inactive_lines, [3 4])','r');
%         plot(lines(active_lines, [1 2])', lines(active_lines, [3 4])','r');
%     end
%     
    lines_votes_vp1 = [lines_votes_vp1(active_lines); lines_votes_vp1(inactive_lines)];
    
    vp(1:2) = computeVPGivenParallelLines(lines(inactive_lines,:));
    
    if DISPLAY_INFO    
        fprintf(fid, 'vp selected = (%.2f, %.2f), vp final = (%.2f, %.2f)\n', vp_selected, vp);
    end
    
    if SHOW_FIGS_PREPROCESS
        figure(3);
       
        plot(vp_selected(1),vp_selected(2),'r*');
        plot(vp(1),vp(2),'b*');        
        plot(lines(inactive_lines, [1 2])', lines(inactive_lines, [3 4])','r');
    end
    vp(1:2) = vp_selected;
end



function vp = computeVPGivenParallelLines(lines)
    num_lines = size(lines,1);
   
    x1 = lines(:, 1);
    x2 = lines(:, 2);
    y1 = lines(:, 3);
    y2 = lines(:, 4);
    
    b = (x2 - x1).*y1 + (y1 - y2).*x1;
    A = zeros(num_lines, 2);
    A(:,1) = (y1-y2);
    A(:,2) = (x2-x1);
    
    vp = A\b;
end