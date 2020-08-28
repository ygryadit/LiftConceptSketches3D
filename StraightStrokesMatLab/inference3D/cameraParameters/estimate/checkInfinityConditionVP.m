function mask_infinite_vps = checkInfinityConditionVP(vp,width)
    w = width;
    h = width;

    mask_infinite_vps = vp(:,1)>10*w | vp(:,2)>10*h;
%     mask_infinite_vps = vp(:,1)>50*w | vp(:,2)>50*h;
end

