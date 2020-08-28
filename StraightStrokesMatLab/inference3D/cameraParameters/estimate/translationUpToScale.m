function t = translationUpToScale(lambda, principal_point_vp, ...
                                      point2DConterpartOrigin, ...
                                      f)
   t(1) = lambda*(point2DConterpartOrigin(1) - principal_point_vp(1))/f;                 
   t(2) = lambda*(point2DConterpartOrigin(2) - principal_point_vp(2))/f;  
   t(3) = lambda;                 
end