function [PosPop, VelPop] = PSO_InitPop(PopNum, ProbDim, DimMinMax)
    PosPop = zeros(PopNum, ProbDim);
    VelPop = zeros(PopNum, ProbDim);
    % Get initial population
    for RowNum = 1:PopNum
        TempPos = zeros(1, ProbDim);
        
        % Randomize Position
        for i = 1:ProbDim
           TempPos(1, i) = (DimMinMax(i, 2) - DimMinMax(i, 1))*rand(1,'double') + DimMinMax(i, 1);
        end
        
        PosPop(RowNum, :) = TempPos; % Add to positions
        
        % Randomize Velocit
        TempVel = zeros(1, ProbDim);
        for i = 1:ProbDim
           VelBound = (DimMinMax(i, 2) - DimMinMax(i, 1))/2;
           TempVel(1, i) = (VelBound - -VelBound)*rand(1,'double') + VelBound;
        end
        VelPop(RowNum, :) = TempVel; % Add to velocities
    end
end

