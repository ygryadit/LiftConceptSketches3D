function cost_full = costPairwise(pCoverage, ...
                                  pIntersections, ...
                                  pAngle, DISP_INFO)

    cost_full = (pCoverage*pIntersections + ...
                  pCoverage*pAngle + ...
                  pAngle*pIntersections + ....
                  pAngle*pIntersections*pCoverage)/4.0;    

% cost_full = (pCoverage*pIntersections + ...
%                   pCoverage*pAngle + ...
%                   pAngle*pIntersections*pCoverage)/3.0;    


    if ~exist('DISP_INFO', 'var')
        DISP_INFO = false;
    end
        
    if DISP_INFO
        fprintf('coverage = %.2f,snap = %.2f,angle = %.2f\n',...
                pCoverage,...
                pIntersections,...
                pAngle);

        fprintf('pAngle*pIntersections = %.2f\n',...
                pAngle*pIntersections);
        fprintf('pCoverage*pIntersections = %.2f\n',...
                pCoverage*pIntersections);
        fprintf('pCoverage*pAngle = %.2f\n',...            
                pCoverage*pAngle);
        fprintf('pAngle*pIntersections*pCoverage = %.2f\n',...            
                pAngle*pIntersections*pCoverage); 

        fprintf('full = %.2f\n', cost_full);
        disp('--------------');
    end
end


%           cost_full =  jointProbablity(...
%                             jointProbablity(...
%                                 jointProbablity(pAngle*pIntersections,...
%                                                 pCoverage*pIntersections),...
%                                 pCoverage*pAngle),...
%                             pAngle*pIntersections*pCoverage);


%         cost_full =  jointProbablity(...
%                             jointProbablity(...
%                                 pCoverage*pIntersections,...
%                                 pCoverage*pAngle),...

%         cost_full = jointProbablity(...
%                                 pCoverage*pIntersections,...
%                                 pCoverage*pAngle);

%         cost_full = jointProbablity(jointProbablity(pCoverage,pAngle),pIntersections); 
%         cost_full =  jointProbablity(...
%                             jointProbablity(...
%                                 jointProbablity(pAngle*pIntersections,...
%                                                 pCoverage*pIntersections),...
%                                 pCoverage*pAngle),...
%                             pAngle*pIntersections*pCoverage);