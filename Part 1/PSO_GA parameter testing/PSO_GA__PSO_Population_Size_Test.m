%LICNACHAN, LANCE OLIVER C.
%2014-64880

format long;

for trials = 1:10
Y = sprintf('Trial: %d',trials);
disp(Y);
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
PopNum = 10;
PSO_Curr=1;
PSO_Max = 1000;
c1 = 1.5;
c2 = 1.5;
wmax = 0.9;
wmin = 0.4;

GA_cross = 0.85;
GA_mut = 0.02;
GA_y = 10;
GA_B = 15;
GA_NumMax = 20;
GA_NumMin = 10;
GA_MinPS = 10;
GA_MaxPS = 25;
GA_MinItr = 50;
GA_MaxItr = 100;
GA_MaxItr = floor(GA_MinItr + ((PSO_Curr/PSO_Max)^GA_B)*(GA_MaxItr-GA_MinItr));
GA_Num = floor(GA_NumMax - ((PSO_Curr/PSO_Max)^GA_y)*(GA_NumMax-GA_NumMin));
GA_PS = floor(GA_MinPS + ((PSO_Curr/PSO_Max)^GA_y)*(GA_MaxPS-GA_MinPS));

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
    
    % Evaluate
    TransFitVal = PSO_GetFitValues(PopNum, TransPos, ProbDim);
    
    % GA Portion
    PSO_Arranged = sort(TransFitVal);
    GA_Num_Curr = 1;
    while GA_Num_Curr <= GA_Num
        % Get one from best individuals
        for RowNum = 1:PopNum
            if TransFitVal(RowNum) == PSO_Arranged(GA_Num_Curr);
               Sel_Indiv = TransPos(RowNum, :); 
            end
        end
        
        % Generate a population with the first indiv being the selected
        % chromosome
        GA_Chroms = GA_InitPop(GA_PS, ProbDim, DimMinMax);
        GA_Chroms(1, :) = Sel_Indiv;
        
        GA_Curr = 1;
        while GA_Curr <= GA_MaxItr
            % Get Fitness
            GA_FitVal = GA_GetFitValues(GA_PS, GA_Chroms, ProbDim);
            TransPop = zeros(GA_PS, ProbDim);
            
            % Keep Elite, Elitist for 50%
            Arranged = sort(GA_FitVal);
            for RowNum = 1:ceil(GA_PS/2)
                for i = 1:GA_PS
                    if Arranged(RowNum) == GA_FitVal(i)
                        TransPop(RowNum, :) = GA_Chroms(i, :);
                    end
                end
            end
            
            % Create Wheel
            GA_RouWheel = GA_CreateWheel(GA_PS, GA_FitVal);
            
            % Create the rest of the population
            for i = ceil(GA_PS/2)+1:GA_PS
                % Select 2 Parents
                [Parent1, Parent2] = GA_Selection(GA_PS, GA_Chroms, ProbDim, GA_FitVal, GA_RouWheel);
                % Cross-over
                SibRep = GA_CrossOver(Parent1, Parent2, GA_cross, ProbDim);
                % Mutate
                if rand() <= GA_mut
                    SibRep = GA_Mutation(SibRep, DimMinMax, ProbDim);
                end
                % Place
                TransPop(i, :) = SibRep;
            end
            
            % Get the best indiv generated
            if GA_Curr == GA_MaxItr
                Best_GA = zeros(1, ProbDim);
                for RowNum = 1:GA_PS
                    if min(GA_FitVal) == GA_FitVal(RowNum)
                        Best_GA(1, :) = GA_Chroms(RowNum, :);
                    end
                end
            end
            
            GA_Chroms = TransPop;
            GA_Curr = GA_Curr + 1;
        end
        % Replace worst individuals
        for RowNum = 1:PopNum
            if TransFitVal(RowNum) == PSO_Arranged((PopNum+1)-GA_Num_Curr);
               TransPos(RowNum,:) = Best_GA(1,:);
            end
        end
        GA_Num_Curr = GA_Num_Curr + 1;
    end
    
    % Update GA_Vars
    GA_MaxItr = floor(GA_MinItr + ((PSO_Curr/PSO_Max)^GA_B)*(GA_MaxItr-GA_MinItr));
    GA_Num = floor(GA_NumMax - ((PSO_Curr/PSO_Max)^GA_y)*(GA_NumMax-GA_NumMin));
    GA_PS = floor(GA_MinPS + ((PSO_Curr/PSO_Max)^GA_y)*(GA_MaxPS-GA_MinPS));
    
    if(PSO_Curr == 1)
        PrevDiff = max(FitVal) - min(FitVal);
    else
        CurrDiff = max(FitVal) - min(FitVal);
        % Check for population convergence
        if PrevDiff - CurrDiff < AcceptThreshold && CurrDiff < AcceptThreshold
            X = sprintf('Population Converged!\nNumber of Iterations: %d\nBest Value: %d\nBest Position: %d %d %d %d %d\n',PSO_Curr,Gbest(ProbDim+1),Gbest(1:ProbDim));
            disp(X);
            break;
        end
        PrevDiff = CurrDiff;
    end
    PosPop = TransPos;
    VelPop = TransVel;
    PSO_Curr = PSO_Curr + 1;
end

if PSO_Curr > PSO_Max
    X = sprintf('Did Not Converge!\nNumber of Iterations: %d\nBest Value: %d\nBest Position: %d %d %d %d %d\n',PSO_Curr,Gbest(ProbDim+1),Gbest(1:ProbDim));
    disp(X);
end

end