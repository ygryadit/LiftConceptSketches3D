% function [failed,vp,p,VoteArrTemp,All_lines] = ...
%             getVPSketch(folder_designer, designer, object_name, view,...
%                        lines_subset_type)
% 
% Output:
%   failed:
%       flag to indicate that the method did not suceed to fing three
%       orthogonal points
%   vp:
%      coordinates of three vanishing points
%   VoteArrTemp:
%       votes of each of the stroke towards one of three vanishing points
%   p:
%       num_lines x 4 (3 vanishign points and the rest)
%       votes normilised by the length of each stroks
% 
%   All_lines:
%       all estiamted strokes that are lines
function [failed,vp,p,...
            vps_selected,... % additional vanishing points for sets of paralel lines
            inds_axis,...    % axis that is orthogonal to the line
            inds_active_lines] = ... % indices of lines that belond to a certain cluster
            getVPSketch(All_lines)
       
    %% Intialise
    global SHOW_FIGS_PREPROCESS;
    global USE_ADDITIONAL_VPs;
    global fid;
    global folderVP;
    failed = false;
    vp = [];
    p = [];
    VoteArrTemp = [];
    
    global sketch_height;
    w = sketch_height; 
    h = sketch_height;
    %% Set display options    
    if SHOW_FIGS_PREPROCESS
        global filepath_sketch_img;
        %img = readImgSketchNew(folder_designer, designer, object_name, view);        
        img = readSketchImg(filepath_sketch_img);
        figure(1);
        imshow(img);
        hold on;

        figure(3);
        plot(0,0);
        hold on;
        imshow(img);
        hold on;
    else
        img =[];
    end
    

    
    %% Convert strokes to a correct represntation:

    All_lines = [All_lines ...
                NaN*ones(size(All_lines(:,1)))...
                NaN*ones(size(All_lines(:,1)))...
                sqrt(((All_lines(:,1)-All_lines(:,2)).^2+(All_lines(:,3)-All_lines(:,4)).^2))];

    
    if (size(All_lines, 1) < 3)
        failed = true;
        return;
    end
     
    %% Compute first vanishing point:
    fprintf(fid, 'II. Computing first vp...');
    [vp, lines_votes_vp1, active_lines, inactive_lines] = getDominantVP(All_lines);
    lines = All_lines(active_lines,:);
    fprintf(fid, ' done.\n');
   
      
%     figure(3);
%     subplot(1,2,1);
%     plot(vp(1), vp(2), '*');
%     hold on;
%     imshow(img);
%     hold on;
%     for i = 1:length(lines_votes_vp1)
%        plot(All_lines(i, 1:2),All_lines(i, 3:4), 'Color', [(1-lines_votes_vp1(i)) 0.0 0.0]);  
%     end
%     
%     
%     subplot(1,2,2);
%     plot(vp(1), vp(2), '*');
%     hold on;
%     imshow(img);
%     hold on;
%     for i = 1:length(lines_votes_vp1)
%        plot(All_lines(i, 1:2),All_lines(i, 3:4), 'Color', [(1-lines_votes_vp1(i)*max_length./All_lines(i,7)) 0.0 0.0]);  
%     end
%     
   
   
    
    %% Work with the remaining lines
    fprintf(fid, 'III. Process remaining lines...\n');
    Xpnts = computeIntersectionPoints(lines);
    inds = find(~isnan(Xpnts(:,1)) & ~isnan(Xpnts(:,2)) & ...
        ~isinf(Xpnts(:,1)) & ~isinf(Xpnts(:,2)));
    Xpnts = Xpnts(inds,:);
    
    if SHOW_FIGS_PREPROCESS
        figure(4);
        plot(0,0);
        hold on;
        imagesc(img);
        axis ij;
        hold on;
        plot(Xpnts(:,1), Xpnts(:,2), 'b*');
        axis ij;
    end
    
     %% Remove some of the points
    fprintf(fid, '\t III.I Clean candidate vps...');
    [Xpnts] = removeRedundantPoints(Xpnts,w,h);
    fprintf(fid, '\t done.\n');
    
    if SHOW_FIGS_PREPROCESS
        figure(4);
        plot(Xpnts(:,1), Xpnts(:,2), 'ro');
    end
    
    %%
    VoteArr = computeLinesPointsVotes([lines;All_lines(inactive_lines,:)],Xpnts);
    Vote=sum(VoteArr(1:size(lines,1),:),1);
  
    
    %Sort vanishign point and their votes:
    [vv, ii]=sort(Vote,'descend');
    Vote = vv(:);
    Xpnts=Xpnts(ii,:);
    VoteArr = VoteArr(:,ii);
    
   
    
    %% Vectorized orthogonality check
    fprintf(fid, '\t III.II Orthogonality check...');
    [pts2,pts1]=find(~triu(ones(length(Vote))));
    npts=length(pts1);
    orthochks=[];
    for pt=1:100000:npts
        tempinds = [pt:min(pt+100000-1,npts)];
        temp_orthochks = checkOthrogonality(ones(length(tempinds),1)*vp(1:2),...
                                            Xpnts(pts1(tempinds),:),...
                                            Xpnts(pts2(tempinds),:),...
                                            w,h);
        orthochks = [orthochks;temp_orthochks(:)];
    end
    orthos = find(orthochks);
    pts1 = pts1(orthos);
    pts2 = pts2(orthos);
    npts = length(pts1);
    fprintf(fid, '\t done.\n');
    
    %% Total vote computation for these points
    totVote = zeros(npts,1);
    for ln=1:length(lines_votes_vp1)
        Votes = [lines_votes_vp1(ln)*ones(npts,1) VoteArr(ln,pts1)' VoteArr(ln,pts2)'];
        Votes = max(Votes,[],2); % Only one line votes for the vanishing point, the one with maximum score
        totVote = totVote+Votes;
    end
    totVote = [pts1(:) pts2(:) totVote(:)];
%     lines = All_lines;

    if size(totVote,1) > 0
        max_length = max(All_lines(:,7));
        [vv ii]=sort(totVote(:,3),'descend');
        vp(3:4) = Xpnts(totVote(ii(1),1),:);
        vp(5:6) = Xpnts(totVote(ii(1),2),:);

        VoteArrTemp = computeLinesPointsVotes(All_lines,[vp(1) vp(2);vp(3) vp(4);vp(5) vp(6)]);
        p=[VoteArrTemp.*max_length./repmat(All_lines(:,7),[1 3]) zeros(size(All_lines,1),1)];%4th vp is outliers
        ind=find(max(p(:,1:3),[],2)< 0.5);
        p(ind,4)=1;
        p=p./repmat(sum(p,2),[1 4]);
%           
%         [~, linemem] = max(p,[],2);
%         
%         display_vp(do_visualize,filepath_vp,img,vp,p,VoteArrTemp,All_lines);
%         
%         grp2=find(linemem==2);
%         grp3=find(linemem==3);
%         vp(3:4) = computeVPGivenParallelLines(All_lines(grp2, :));
%         vp(5:6) = computeVPGivenParallelLines(All_lines(grp3, :));
%         
        %     [vv linemem] = max(VoteArrTemp,[],2);
    
        
        % Plot three vanishing points:
        if SHOW_FIGS_PREPROCESS
            visualizeVP(SHOW_FIGS_PREPROCESS,folderVP,img,vp,p,All_lines);
        end
    else
        failed = true;
    end
    
    %% Additional sets of lines with priors:
    vps_selected = [];
    inds_axis = [];
    inds_active_lines = [];
    
    if failed        
        return;
    end
    if USE_ADDITIONAL_VPs
    [vps_selected,...
     inds_axis,...
     inds_active_lines] = computeRotatedSetsOfParallelLines(vp,p,All_lines,w,h,img);
    end
end
