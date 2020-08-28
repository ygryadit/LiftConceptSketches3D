% function [orthochk]=chckothrogonalityvector(vp1s,vp2s,vp3s,w,h,center)
 function [orthochk]=checkOthrogonality(vp1s,vp2s,vp3s,w,h)
% vp1s, vp2s, and vp3s multiple vanishing points.

orthochk=zeros(size(vp1s,1),1);

%inf conditions
inf1 = vp1s(:,1)>10*w | vp1s(:,2)>10*h;
inf2 = vp2s(:,1)>10*w | vp2s(:,2)>10*h;
inf3 = vp3s(:,1)>10*w | vp3s(:,2)>10*h;

inds_fff = find(~inf1 & ~inf2 & ~inf3);
% Rearrange to have third point at in infinity
inds = find(inf1 & ~inf2 & ~inf3);
temp = vp1s(inds,:);
vp1s(inds,:) = vp3s(inds,:);
vp3s(inds,:) = temp;
inds = find(~inf1 & inf2 & ~inf3);
temp = vp2s(inds,:);
vp2s(inds,:) = vp3s(inds,:);
vp3s(inds,:) = temp;
% We moved infinite point to the third row
inds_ffi = find((inf1+inf2+inf3)==1); % Only one point is at infinity


inds = find(inf1 & ~inf2 & inf3);
% Rearrange to make 2nd and 3rd point being at infinity
temp = vp2s(inds,:);
vp2s(inds,:) = vp1s(inds,:);
vp1s(inds,:) = temp;
inds = find(inf1 & inf2 & ~inf3);
temp = vp3s(inds,:);
vp3s(inds,:) = vp1s(inds,:);
vp1s(inds,:) = temp;
inds_fii = find((inf1+inf2+inf3)==2);
% inds_fii = find(~inf1 & inf2 & inf3);
inds_iii = find(inf1 & inf2 & inf3);

if numel(inds_fff)>0
    fff_orthochk = orthochk_3finite_vp(vp1s,vp2s,vp3s,inds_fff,w,h);
    orthochk(inds_fff) = fff_orthochk;
end

if numel(inds_ffi)>0
%     v1sf = vp1s(inds_ffi,:);
%     v2sf = vp2s(inds_ffi,:);
%     v3si = vp3s(inds_ffi,:);
% 
%     temp_orthochk = false(size(inds_ffi));
%     
%     % In the original Rother's paper the point is selected lying on the
%     % line between thw two finite vanishing points and closer to the image
%     % center. 
%     % In our case it makes sense to allign the center with the cneter mass
%     % of the sketched shape itself.
%     
%     
% %     r=((w/2-v1sf(:,1)).*(v2sf(:,1)-v1sf(:,1))+(h/2-v1sf(:,2)).*(v2sf(:,2)-v1sf(:,2)))./((v2sf(:,1)-v1sf(:,1)).^2+(v2sf(:,2)-v1sf(:,2)).^2);    
% 
%     r=((center(:,1)-v1sf(:,1)).*(v2sf(:,1)-v1sf(:,1))+(center(:,2)-v1sf(:,2)).*(v2sf(:,2)-v1sf(:,2)))./((v2sf(:,1)-v1sf(:,1)).^2+(v2sf(:,2)-v1sf(:,2)).^2);    
%     
%     u0= v1sf(:,1) + r.*(v2sf(:,1)-v1sf(:,1));
%     v0= v1sf(:,2) + r.*(v2sf(:,2)-v1sf(:,2));
%     ind_feasible_principal_point = r > 0 & r < 1 & u0 <= 0.7*w & u0 >= 0.3*w & v0 <= 0.7*h & v0 >= 0.3*h;
%     
% %     %Line passign through vanishing points:
% %     slope = (v2sf(:,2)-v1sf(:,2))./(v2sf(:,1)-v1sf(:,1));
% %     
% %     %y = slope*x + b1
% %     b1 =  v1sf(:,2) - v1sf(:,1).*slope;
% %     
% %     %y = -1/slope*x + b2  -- orthogonal line passing through the center of
% %     %the drawing:
% %     b2 =  center(:,2) + center(:,1)./slope;
% %     
% %     % Obtain center as the intersection point
% %     u0 = (b2 - b1)./(slope + 1.0./slope);
% %     v0 = slope.*u0 + b1;
%     
% %     
%     
% %     inds = find(u0 > 0.7*w | u0 < 0.3*w | v0 > 0.7*h | v0 < 0.3*h | isnan(u0) | isnan(v0));
% %     temp_orthochk(inds) = 0;
% 
%     vec1 = v1sf - principal_point;
%     vec2 = principal_point - v2sf;    
%     f = sqrt(abs(dot(vec1,vec2,2))); 
%     fov = 2*atand(w./(2.0*f));
%     inf_feasible_fov = fov > 10 & fov < 120; 
% 
% 
%     temp=u0.*(v1sf(:,1)+v2sf(:,1))+v0.*(v2sf(:,2)+v1sf(:,2))-(v1sf(:,1).*v2sf(:,1)+v2sf(:,2).*v1sf(:,2)+u0.^2+v0.^2);
% %     inds = find(temp<0 | isnan(temp));
% %     temp_orthochk(inds) = 0;
%     f = (temp).^(0.5);
% %     inds = find(f<=0 | f>=5000 | isnan(f));
% %     temp_orthochk(inds) = 0;
% 
%     vec1=[v2sf(:,1)-v1sf(:,1) v2sf(:,2)-v1sf(:,2)]; 
%     vec2=[v3si(:,1) v3si(:,2)];
%     dot12 = sum(vec1.*vec2,2);
%     norm1 = (sum(vec1.*vec1,2).^.5);
%     norm2 = (sum(vec2.*vec2,2).^.5);
% %     inds = find(abs(dot12./norm1./norm2)>=0.1 | isnan(dot12) | isnan(norm1) | isnan(norm2));
% %     temp_orthochk(inds) = 0;
% 
% 
%     inds=find(r > 0 & r < 1 & u0 <= 0.7*w & u0 >= 0.3*w & v0 <= 0.7*h & v0 >= 0.3*h & temp > 0 & f> 0 & f<=5000 ...
%      & abs(dot12./norm1./norm2)< 0.1);
%     temp_orthochk(inds) = true;

    %     temp_orthochk2 = zeros(size(temp_orthochk));
    %     for i=1:length(temp_orthochk2)
    %         temp_orthochk2(i) = chckothrogonality([v1sf(i,:) v2sf(i,:) v3sf(i,:)],w,h);
    %     end
    ffi_orthochk = orthochk_2finite_vp(vp1s,vp2s,vp3s,inds_ffi,w,h);
    orthochk(inds_ffi) = ffi_orthochk;
end

if numel(inds_fii)>0
    v1sf = vp1s(inds_fii,:);
    v2si = vp2s(inds_fii,:);
    v3si = vp3s(inds_fii,:);
%     temp_orthochk = ones(size(inds_fii));

    temp_orthochk = zeros(size(inds_fii));

    vec1 = [v2si(:,1) v2si(:,2)];
    vec2 = [v3si(:,1) v3si(:,2)];
    dot12 = sum(vec1.*vec2,2);
    norm1 = (sum(vec1.*vec1,2).^.5);
    norm2 = (sum(vec2.*vec2,2).^.5);
    
%     inds = find(abs(dot12./norm1./norm2)>=0.1 | isnan(dot12) | isnan(norm1) | isnan(norm2));
%     temp_orthochk(inds) = 0;
% 
%     inds = find(v1sf(:,1)<0.3*w | v1sf(:,1)>0.7*w | v1sf(:,2)<0.3*h | v1sf(:,2)>0.7*h | isnan(v1sf(:,1)) | isnan(v1sf(:,2)));
%     temp_orthochk(inds) = 0;

%     inds = find(abs(dot12./norm1./norm2)<=0.1 & v1sf(:,1)> 0 & v1sf(:,1)<=0.7*w & v1sf(:,2)> 0 & v1sf(:,2)<=0.7*h );
    inds = abs(dot12./norm1./norm2)<=0.1 & ind_feasible_principle_point_mild(v1sf(:,1), v1sf(:,2), w, h) ;
%     inds_orthogonal = checkOrthogonalityCriterion(v1sf, v2sf, v3sf);
    
%     inds = inds_orthogonal & inds;

    temp_orthochk(inds) = true;


    %     temp_orthochk2 = zeros(size(temp_orthochk));
    %     for i=1:length(temp_orthochk2)
    %         temp_orthochk2(i) = chckothrogonality([v1sf(i,:) v2sf(i,:) v3sf(i,:)],w,h);
    %     end

    orthochk(inds_fii) = temp_orthochk;
end

if numel(inds_iii)>0
    orthochk(inds_iii) = 0;
end

end

function inds = checkOrthogonalityCriterion(v1sf, v2sf, v3sf)
    side1 = (v1sf - v2sf);
    side2 = (v1sf - v3sf);
    side3 = (v2sf - v3sf);
    
    side1 = side1./normArray(side1);
    side2 = side2./normArray(side2);
    side3 = side3./normArray(side3);
    
    theta1 = acos(dot(side1,side2,2));
    theta2 = acos(dot(-side3,-side2,2));
    theta3 = acos(dot(-side1,side3,2));

    inds = theta1 < pi/2 & theta2 < pi/2 & theta3 < pi/2;
end


function n = normArray(array)
% Input: 
%   array: nx2 array
% 
% Output:
%   n: the norm along first dimension
   
    n = (sqrt(sum((array').^2)))';
end

function principal_point = findPrincipalPointFocalLength3VP(v1sf, v2sf, v3sf)
% Find principal point given triplets of three finite vanishing points.
    Mats_11 = v1sf(:,1)+v2sf(:,1);
    Mats_12 = v1sf(:,2)+v2sf(:,2);
    Mats_13 = v1sf(:,1).*v2sf(:,1)+v1sf(:,2).*v2sf(:,2);
    Mats_21 = v1sf(:,1)+v3sf(:,1);
    Mats_22 = v1sf(:,2)+v3sf(:,2);
    Mats_23 = v1sf(:,1).*v3sf(:,1)+v1sf(:,2).*v3sf(:,2);
    Mats_31 = v3sf(:,1)+v2sf(:,1);
    Mats_32 = v3sf(:,2)+v2sf(:,2);
    Mats_33 = v3sf(:,1).*v2sf(:,1)+v3sf(:,2).*v2sf(:,2);

    A_11 = Mats_11-Mats_21; A_12 = Mats_12-Mats_22;
    A_21 = Mats_11-Mats_31; A_22 = Mats_12-Mats_32;
    b_1 = Mats_13-Mats_23; b_2 = Mats_13-Mats_33;
    detA = A_11.*A_22-A_12.*A_21;
    principal_point = zeros(length(v1sf), 2);
    principal_point(:,1) = (A_22.*b_1-A_12.*b_2)./detA;
    principal_point(:,2) = (A_11.*b_2-A_21.*b_1)./detA;
end

function fff_orthochk = orthochk_3finite_vp(vp1s,vp2s,vp3s,inds_fff,w,h)
    fff_orthochk = false(size(inds_fff));    
    
    v1sf = vp1s(inds_fff,:);
    v2sf = vp2s(inds_fff,:);
    v3sf = vp3s(inds_fff,:);
    
    %I. Check orthogonality condition:
%     inds_orthogonal  = 1:length(inds_fff);
    inds_orthogonal = find(checkOrthogonalityCriterion(v1sf, v2sf, v3sf));
    v1sf = v1sf(inds_orthogonal, :);
    v2sf = v2sf(inds_orthogonal, :);
    v3sf = v3sf(inds_orthogonal, :);
        
    %Check the rest of conditions only on those points that passed intial
    %test on orthogonality:
    
    %II. Find principal point:
    principal_point = findPrincipalPointFocalLength3VP(v1sf, v2sf, v3sf);
    u0 = principal_point(:,1);
    v0 = principal_point(:,2);
    ind_feasible_principal_point = ind_feasible_principle_point(u0, v0, w, h); %ind_feasible_principle_point_mild(u0, v0, w, h);
    
    %III. Find focal length and check that it is feasible:
    vec1 = v1sf - principal_point;
    vec2 = principal_point - v2sf;    
    f = sqrt(abs(dot(vec1,vec2,2))); 
    fov = 2*atand(w./(2.0*f));
    inf_feasible_fov = fov > 10 & fov < 120; 
    
%     inds = find(u0 <= 0.7*w & u0 >= 0.3*w & v0 <= 0.7*h & v0 >= 0.3*h & f> 0 & f<=5000 & temp> 0);

 
    inds = ind_feasible_principal_point & inf_feasible_fov;
    inds = inds_orthogonal(inds);    
    fff_orthochk(inds) = true;
end

function ffi_orthochk = orthochk_2finite_vp(vp1s,vp2s,vp3s,inds_ffi,w,h)
    ffi_orthochk = false(size(inds_ffi));

    v1sf = vp1s(inds_ffi,:);
    v2sf = vp2s(inds_ffi,:);
    v3si = vp3s(inds_ffi,:);

    % I. Find principal point and check its feasibility:
    % In the original Rother's paper the point is selected lying on the
    % line between the two finite vanishing points and closer to the image
    % center. 
    % In our case it makes sense to allign the center with the cneter mass
    % of the sketched shape itself.
      
%     r=((w/2-v1sf(:,1)).*(v2sf(:,1)-v1sf(:,1))+(h/2-v1sf(:,2)).*(v2sf(:,2)-v1sf(:,2)))./((v2sf(:,1)-v1sf(:,1)).^2+(v2sf(:,2)-v1sf(:,2)).^2);    
    
%     ref_point = center;
    ref_point = v3si;
    
    % Find projection of the infinite vanishing point on the line
    % connecting the two vanishing points:
    r=((ref_point(:,1)-v1sf(:,1)).*(v2sf(:,1)-v1sf(:,1))+(ref_point(:,2)-v1sf(:,2)).*(v2sf(:,2)-v1sf(:,2)))...
        ./((v2sf(:,1)-v1sf(:,1)).^2+(v2sf(:,2)-v1sf(:,2)).^2);    
    
    u0= v1sf(:,1) + r.*(v2sf(:,1)-v1sf(:,1));
    v0= v1sf(:,2) + r.*(v2sf(:,2)-v1sf(:,2));
    
    point_is_between_two_finite_points = r > 0 & r < 1;
    ind_feasible_principal_point = point_is_between_two_finite_points & ind_feasible_principle_point_2finitevp(u0, w);
    
    principal_point = zeros(size(inds_ffi));
    principal_point(:,1) = u0;
    principal_point(:,2) = v0;
    
    % II. Find focal length and check its feasibility:
    vec1 = v1sf - principal_point;
    vec2 = principal_point - v2sf;    
    f = sqrt(abs(dot(vec1,vec2,2))); 
    fov = 2*atand(w./(2.0*f));
    inf_feasible_fov = fov > 10 & fov < 120; 
    
    inds = ind_feasible_principal_point & inf_feasible_fov;

%     % III. Check oorthognality criteria
%     vec1=[v2sf(:,1)-v1sf(:,1) v2sf(:,2)-v1sf(:,2)]; 
%     vec2=[v3si(:,1) v3si(:,2)];
%     dot12 = sum(vec1.*vec2,2);
%     norm1 = (sum(vec1.*vec1,2).^.5);
%     norm2 = (sum(vec2.*vec2,2).^.5);
% %     inds = find(abs(dot12./norm1./norm2)>=0.1 | isnan(dot12) | isnan(norm1) | isnan(norm2));
% %     ffi_orthochk(inds) = 0;
% 
% 
%     inds= ind_feasible_principal_point & inf_feasible_fov...
%      & abs(dot12./norm1./norm2)< 0.1;
%  
 
    ffi_orthochk(inds) = true;

end

function ind = ind_feasible_principle_point(u0, v0, w, h)
%     ind = u0 <= 3*w & u0 >= -3*w & v0 <= 3*h & v0 >= -3*h;
    ind = u0 <= 0.7*w & u0 >= 0.3*w & v0 <= 0.7*h & v0 >= 0.3*h;
end

function ind = ind_feasible_principle_point_mild(u0, v0, w, h)
%     ind = u0 <= 2*w & u0 >= -2*w & v0 <= 2*h & v0 >= -2*h;
    ind = u0 <= 0.7*w & u0 >= 0.3*w & v0 <= 0.7*h & v0 >= 0.3*h;
end

function ind = ind_feasible_principle_point_2finitevp(u0, w)
    %For now only width is considered since in case of vertical infinite
    %vanishing point the h can be out of the omage plane for the principle
    %point.
    ind = u0 <= 0.7*w & u0 >= 0.3*w;
end