% Function that produces the initial PSO population
function [PosPop, VelPop] = PSO_InitPop(PopNum, ProbDim, DimMinMax)
    PosPop = zeros(PopNum, ProbDim);
    VelPop = zeros(PopNum, ProbDim);
    
    % Create each particle
    for RowNum = 1:PopNum
        TempPos = zeros(1, ProbDim);
        
        % Randomize Position
        for i = 1:ProbDim
           TempPos(1, i) = (DimMinMax(i, 2) - DimMinMax(i, 1))*rand(1,'double') + DimMinMax(i, 1);
        end
%         uniq=0;
%         % Check for uniqueness
%         if length(TempPos) == length(unique(TempPos))
%             uniq=1;
%         end
%         % Resolve uniqueness
%         TempPos = Resolve_Uniqueness(TempPos, ProbDim, uniq);
        
        % Add to positions
        PosPop(RowNum, :) = TempPos; 
        
        % Randomize Velocity
        TempVel = zeros(1, ProbDim);
        for i = 1:ProbDim
           % Set Boundary
           VelBound = (DimMinMax(i, 2) - DimMinMax(i, 1));
           % Generate within bounds
           TempVel(1, i) = (VelBound - -VelBound)*rand(1,'double') - VelBound;
        end 
        % Add to velocities
        VelPop(RowNum, :) = TempVel;
    end
end

