function [Par, is_singular] = CircleFitByPratt(XY)

%--------------------------------------------------------------------------
%  
%     Algebraic circle fit by Pratt
%      V. Pratt, "Direct least-squares fitting of algebraic surfaces",
%      Computer Graphics, Vol. 21, pages 145-152 (1987)
%
%     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
%
%     Output: Par = [a b R] is the fitting circle:
%                           center (a,b) and radius R
%
%     Note: this is a version optimized for stability, not for speed
%
%--------------------------------------------------------------------------

centroid = mean(XY);   % the centroid of the data set

[U,S,V]=svd([(XY(:,1)-centroid(1)).^2+(XY(:,2)-centroid(2)).^2,...
        XY(:,1)-centroid(1), XY(:,2)-centroid(2), ones(size(XY,1),1)],0);
is_singular = false;
if (S(4,4)/S(1,1) < 1e-12)   %  singular case
    A = V(:,4);
%     disp('Pratt singular case');
    is_singular = true;
else                         %  regular case
    W=V*S;
    Binv = [0 0 0 -0.5; 0 1 0 0; 0 0 1 0; -0.5 0 0 0];
    [E,D] = eig(W'*Binv*W);
    [Dsort,ID] = sort(diag(D));
    A = E(:,ID(2));
    for i=1:4
        S(i,i)=1/S(i,i); 
    end
    A = V*S*A;
end
% tt = A(1)*(XY(:,1)-centroid(1)).^2+(XY(:,2)-centroid(2)).^2 +...
% A(2)*(XY(:,1)-centroid(1)) + A(3)*(XY(:,2)-centroid(2)) + A(4)*ones(size(XY,1),1);  


Par = -(A(2:3))'/A(1)/2 + centroid;
Par = [Par , sqrt(A(2)^2+A(3)^2-4*A(1)*A(4))/abs(A(1))/2];


tt = sqrt((XY(:,1) - Par(1)).^2 +  (XY(:,2) - Par(2)).^2) - Par(3);
Par(4) = sum(abs(tt))./length(tt);
if is_singular
    Par(4) = 0;
end

end   %  PrattSVD



% % https://people.cas.uab.edu/~mosya/cl/MATLABcircle.html
% %--------------------------------------------------------------------------
% %  
% %     Circle fit by Pratt
% %      V. Pratt, "Direct least-squares fitting of algebraic surfaces",
% %      Computer Graphics, Vol. 21, pages 145-152 (1987)
% %
% %     Input:  XY(n,2) is the array of coordinates of n points x(i)=XY(i,1), y(i)=XY(i,2)
% %
% %     Output: Par = [a b R] is the fitting circle:
% %                           center (a,b) and radius R
% %
% %     Note: this fit does not use built-in matrix functions (except "mean"),
% %           so it can be easily programmed in any programming language
% %
% %--------------------------------------------------------------------------
% exit_flag =0;
% n = size(XY,1);      % number of data points
% centroid = mean(XY);   % the centroid of the data set
% %     computing moments (note: all moments will be normed, i.e. divided by n)
% Mxx=0; Myy=0; Mxy=0; Mxz=0; Myz=0; Mzz=0;
% for i=1:n
%     Xi = XY(i,1) - centroid(1);  %  centering data
%     Yi = XY(i,2) - centroid(2);  %  centering data
%     Zi = Xi*Xi + Yi*Yi;
%     Mxy = Mxy + Xi*Yi;
%     Mxx = Mxx + Xi*Xi;
%     Myy = Myy + Yi*Yi;
%     Mxz = Mxz + Xi*Zi;
%     Myz = Myz + Yi*Zi;
%     Mzz = Mzz + Zi*Zi;
% end
%    
% Mxx = Mxx/n;
% Myy = Myy/n;
% Mxy = Mxy/n;
% Mxz = Mxz/n;
% Myz = Myz/n;
% Mzz = Mzz/n;
% %    computing the coefficients of the characteristic polynomial
% Mz = Mxx + Myy;
% Cov_xy = Mxx*Myy - Mxy*Mxy;
% Mxz2 = Mxz*Mxz;
% Myz2 = Myz*Myz;
% A2 = 4*Cov_xy - 3*Mz*Mz - Mzz;
% A1 = Mzz*Mz + 4*Cov_xy*Mz - Mxz2 - Myz2 - Mz*Mz*Mz;
% A0 = Mxz2*Myy + Myz2*Mxx - Mzz*Cov_xy - 2*Mxz*Myz*Mxy + Mz*Mz*Cov_xy;
% A22 = A2 + A2;
% epsilon=1e-12; 
% ynew=1e+20;
% IterMax=20;
% xnew = 0;
% %    Newton's method starting at x=0
% for iter=1:IterMax
%     yold = ynew;
%     ynew = A0 + xnew*(A1 + xnew*(A2 + 4.*xnew*xnew));
%     if (abs(ynew)>abs(yold))
%         disp('Newton-Pratt goes wrong direction: |ynew| > |yold|');
%         exit_flag =1;
%         xnew = 0;
%         break;
%     end
%     Dy = A1 + xnew*(A22 + 16*xnew*xnew);
%     xold = xnew;
%     xnew = xold - ynew/Dy;
%     
%     if (abs((xnew-xold)/xnew) < epsilon), break, end
%     if (iter >= IterMax)
%         disp('Newton-Pratt will not converge');
%         exit_flag =2;
%         xnew = 0;
%     end
%     if (xnew<0.)
%         fprintf(1,'Newton-Pratt negative root:  x=%f\n',xnew);
%         exit_flag =3;
%         xnew = 0;
%     end
% end
% 
% %    computing the circle parameters
% DET = xnew*xnew - xnew*Mz + Cov_xy;
% % fprintf('xnew = %.3f, ynew = %.3f, Dy = %.3f\n', xnew, ynew, Dy);
% Center = [Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy]/DET/2;
% Par = [Center+centroid , sqrt(Center*Center'+Mz+2*xnew), xnew];
% end    %    CircleFitByPratt
