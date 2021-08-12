function [PosPop, VelPop] = PSO_InitPop(PopNum, ProbDim, DimMinMax)
    PosPop = zeros(PopNum, ProbDim);
    VelPop = zeros(PopNum, ProbDim);
    
    for RowNum = 1:PopNum
        TempPos = zeros(1, ProbDim);
        
        % Randomize Position
        for i = 1:ProbDim
           TempPos(1, i) = randi([DimMinMax(i, 1), DimMinMax(i, 2)], 1);
        end
        
        PosPop(RowNum, :) = TempPos; % Add to positions
        
        % Randomize Velocit
        TempVel = zeros(1, ProbDim);
        for i = 1:ProbDim
           VelBound = ceil((DimMinMax(i, 2) - DimMinMax(i, 1))/2);
           TempVel(1, i) = randi([-VelBound, VelBound], 1);
        end
        VelPop(RowNum, :) = TempVel; % Add to velocities
    end
end

