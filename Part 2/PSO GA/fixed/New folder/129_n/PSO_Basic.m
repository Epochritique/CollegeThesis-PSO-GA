%LICNACHAN, LANCE OLIVER C.
%2014-64880

format long;

% Variables specific to the problem
ProbDim = 5;
DimMinMax = zeros(ProbDim, 2);
DimMinMax(1, :) = [0 1];
DimMinMax(1, :) = [0 1];
DimMinMax(1, :) = [0 1];
DimMinMax(1, :) = [0 1];
DimMinMax(1, :) = [0 1];

% Variables specific to the algorithm
AcceptThreshold = 1e-6;
PopNum = 50;
PSO_Curr = 1;
PSO_Max = 1000;
c1 = 1.5;
c2 = 1.5;
wmax = 0.9;
wmin = 0.4;

TransPos = zeros(PopNum, ProbDim);
TransVel = zeros(PopNum, ProbDim);

% Initialization Step
[PosPop, VelPop] = PSO_InitPop(PopNum, ProbDim, DimMinMax);
FitVal = PSO_GetFitValues(PopNum, PosPop, ProbDim);
[Pbest, Gbest]= PSO_GetInitPGBest(PopNum, PosPop, ProbDim, FitVal);

PrevDiff = 0;
while PSO_Curr <= PSO_Max
    % Evaluate
    FitVal = PSO_GetFitValues(PopNum, PosPop, ProbDim);
    
    % Get best values
    [Pbest, Gbest] = PSO_GetPGBest(PopNum, PosPop, ProbDim, FitVal, Pbest, Gbest);
    
    % Change value according to how current iteration
    w = wmax-(wmax-wmin)*(PSO_Curr/PSO_Max);
    
    % Calculate new velocities and move
    [TransPos, TransVel] = PSO_ChangeVel(PopNum, PosPop, VelPop, ProbDim, Pbest, Gbest, w, c1, c2);
    
    if(PSO_Curr == 1)
        PrevDiff = max(FitVal) - min(FitVal);
    else
        CurrDiff = max(FitVal) - min(FitVal);
        % Check for population convergence
        if PrevDiff - CurrDiff < AcceptThreshold && CurrDiff < AcceptThreshold
            X = sprintf('Population Converged!\nNumber of Iterations: %d\nBest Value: %d\nBest Position: %d %d %d %d %d',PSO_Curr,Gbest(ProbDim+1),Gbest(1:ProbDim));
            disp(X);
            break;
        end
        PrevDiff = CurrDiff;
    end
    disp(X);
    PosPop = TransPos;
    VelPop = TransVel;
    PSO_Curr = PSO_Curr + 1;
end

if PSO_Curr > PSO_Max
    X = sprintf('Did Not Converge!\nNumber of Iterations: %d\nBest Value: %d\n',PSO_Curr,Gbest(ProbDim+1));
    disp(X);
end

