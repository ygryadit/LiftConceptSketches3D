%  Computes voting weights as in Varsha Hedau, Derek Hoiem, David Forsyth, 
%  “Recovering the Spatial Layout of Cluttered Rooms,” in the Twelfth IEEE 
%  International Conference on Computer Vision, 2009. Equation 1

function VoteArr = computeLinesPointsVotes(lines, Xpnts)
% ta=pi/3;

ta = pi/3; %threshold on line to belong to the  pi/6 failure cases

max_length=max(lines(:,7));
% w1=0.60;
% w2=0.40;

VoteArr=zeros(size(lines,1),size(Xpnts,1));
for line_num = 1:size(lines,1)
    [theta, vp_outside_line_segment]=computeDistanceLineVPs(lines(line_num,:),Xpnts);    
    indd=find(theta < ta & vp_outside_line_segment);
    % Carsten Rother: "A new approach to vanishing point detection in architectural
    % environments" Eq.1 term 1:
    VoteArr(line_num,indd) = (1-theta(indd)/ta);
    % Carsten Rother: "A new approach to vanishing point detection in architectural
    % environments" Eq.1
    %     VoteArr(line_num,indd)= w1*(1-(1/ta).*theta(indd))+ w2*lines(line_num,7)/max_length;
end
VoteArr = exp(-(1-VoteArr).^2/0.1/0.1/2);
VoteArr = VoteArr.*repmat(lines(:,7),1,size(VoteArr,2))/max_length;
return;
