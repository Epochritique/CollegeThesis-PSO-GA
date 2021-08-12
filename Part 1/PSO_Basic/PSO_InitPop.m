function [PosPop, VelPop] = PSO_InitPop(PopNum, ProbDim, DimMinMax)
    PosPop = zeros(PopNum, ProbDim);
    VelPop = zeros(PopNum, ProbDim);
    % Get initial population
    for RowNum = 1:PopNum
        TempPos = zeros(1, ProbDim);
        
        % Randomize the position
        for i = 1:ProbDim
            TempPos(1, i) = randi(2)-1;
            % TempPos(1, i) = (DimMinMax(i,2) - DimMinMax(i,1))*rand(1,'double') + DimMinMax(i,1);
        end
        PosPop(RowNum, :) = TempPos; % Add to positions
        
        % Randomize the velocity
        TempVel = zeros(1, ProbDim);
        for i = 1:ProbDim
            InitVelo = (31/2);
            % InitVelo = ((DimMinMax(i, 2) - DimMinMax(i, 1))/2);
            TempVel(1, i) = (InitVelo - -InitVelo)*rand(1,'double') + InitVelo;
        end
        VelPop(RowNum, :) = TempVel; % Add to velocities
    end
end

