%LICNACHAN, LANCE OLIVER C.
%2014-64880

format long;
rng('shuffle');

ProbDim = 4;
ConsNum = 1;
Ans = zeros(30, ProbDim+ConsNum+1);
for trials = 1:30
Y = sprintf('Trial: %d',trials);
disp(Y);

ctr=0;

% Variables specific to the problem

DimMinMax = zeros(ProbDim, 2);
DimMinMax(1, :) = [12 60];
DimMinMax(2, :) = [12 60];
DimMinMax(3, :) = [12 60];
DimMinMax(4, :) = [12 60];

% Variables specific to the algorithm
AcceptThreshold = 1e-6;
PopNum = 20*ProbDim;
PSO_Curr = 1;
PSO_Max = 500;
c1 = 1.5;
c2 = 1.5;
wmax = 0.9;
wmin = 0.4;

GA_cross = 0.85;
GA_mut = 0.02;
GA_y = 10;
GA_B = 15;
GA_NumMax = 20;
GA_NumMin = 1;
GA_MinPS = 10;
GA_MaxPS = 5;
GA_MinItr = 10;
GA_MaxItr = 15;
GA_MaxItr = floor(GA_MinItr + ((PSO_Curr/PSO_Max)^GA_B)*(GA_MaxItr-GA_MinItr));
GA_Num = floor(GA_NumMax - ((PSO_Curr/PSO_Max)^GA_y)*(GA_NumMax-GA_NumMin));
GA_PS = floor(GA_MinPS + ((PSO_Curr/PSO_Max)^GA_y)*(GA_MaxPS-GA_MinPS));

TransPos = zeros(PopNum, ProbDim);
TransVel = zeros(PopNum, ProbDim);

% Initialization Step
[PosPop, VelPop] = PSO_InitPop(PopNum, ProbDim, DimMinMax);
FitVal = PSO_GetFitValues(PopNum, PosPop, ProbDim, DimMinMax);
[Pbest, Gbest]= PSO_GetInitPGBest(PopNum, PosPop, ProbDim, FitVal);

PrevDiff = 0;
while PSO_Curr <= PSO_Max
%     disp(Gbest(ProbDim+1));
    
    % Evaluate
    FitVal = PSO_GetFitValues(PopNum, PosPop, ProbDim, DimMinMax);
% 
%     clf;    %clear frame
%     figure(1);
%     hold on;
%     posit = 1:PopNum;
%     plot(posit,FitVal,'.r','MarkerSize', 10);
%     M(PSO_Curr)=getframe(gca);
    
    % Get best values
    [Pbest, Gbest] = PSO_GetPGBest(PopNum, PosPop, ProbDim, FitVal, Pbest, Gbest);
    
    if(PSO_Curr == PSO_Max)
        break;
    end
    % Change value according to how current iteration
    w = wmax-(wmax-wmin)*(PSO_Curr/PSO_Max);
    
    % Calculate new velocities and move
    [TransPos, TransVel] = PSO_ChangeVel(PopNum, PosPop, VelPop, ProbDim, Pbest, Gbest, w, c1, c2, DimMinMax);
    
    % Evaluate
    TransFitVal = PSO_GetFitValues(PopNum, TransPos, ProbDim, DimMinMax);
    
    % GA Portion
    PSO_Arranged = sort(TransFitVal);
    GA_Num_Curr = 1;
    while GA_Num_Curr <= GA_Num
        % Get one from best individuals
        for RowNum = 1:PopNum
            if TransFitVal(RowNum) == PSO_Arranged(GA_Num_Curr);
               Sel_Indiv = TransPos(RowNum, :);
               break;
            end
        end
        
        % Generate a population with the first indiv being the selected
        % chromosome
        GA_Chroms = GA_InitPop(GA_PS, ProbDim, DimMinMax);
        GA_Chroms(1, :) = Sel_Indiv;
        
        GA_Fit_Elite = PSO_Arranged(GA_Num_Curr);
        GA_Fit_Chrom = Sel_Indiv;
        GA_Curr = 1;
        while GA_Curr <= GA_MaxItr
            % Get Fitness
            GA_FitVal = PSO_GetFitValues(GA_PS, GA_Chroms, ProbDim, DimMinMax);
            TransPop = zeros(GA_PS, ProbDim);
            
            % Keep Elite
            Arranged = sort(GA_FitVal);
            if Arranged(1) < GA_Fit_Elite
                GA_Fit_Elite = Arranged(1);
                for i = 1:GA_PS
                    if Arranged(1) == GA_FitVal(i)
                        GA_Fit_Chrom = GA_Chroms(i,:);
                    end
                end
            end
            
            % Create Wheel
            GA_RouWheel = GA_CreateWheel(GA_PS, GA_FitVal);
            
            % Create the population
            for i = 1:GA_PS
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
            
            GA_Chroms = TransPop;
            GA_Curr = GA_Curr + 1;
        end
        % Replace the individual
        for RowNum = 1:PopNum
            if TransFitVal(RowNum) == PSO_Arranged(GA_Num_Curr);
                TransPos(RowNum,:) = GA_Fit_Chrom(1,:);
                break;
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
        % disp(CurrDiff);
        % Check for population convergence
        GearRatio = (Gbest(1)*Gbest(2))/(Gbest(3)*Gbest(4));
        if PrevDiff - CurrDiff < AcceptThreshold && CurrDiff < AcceptThreshold
            X = sprintf('Population Converged!\nNumber of Iterations: %0.15f\nBest Value: %0.15f\nBest Position: %0.15f %0.15f %0.15f %0.15f\nGear Ratio: %0.15f\n',PSO_Curr,Gbest(ProbDim+1),Gbest(1:ProbDim),GearRatio);
            disp(X);
            break;
        end
        PrevDiff = CurrDiff;
    end
    
    PosPop = TransPos;
    VelPop = TransVel;
    PSO_Curr = PSO_Curr + 1;
end

if PSO_Curr >= PSO_Max
    GearRatio = (Gbest(1)*Gbest(2))/(Gbest(3)*Gbest(4));
    X = sprintf('Did Not Converge!\nNumber of Iterations: %0.15f\nBest Value: %0.15f\nBest Position: %0.15f %0.15f %0.15f %0.15f\nGear Ratio: %0.15f\n',PSO_Curr,Gbest(ProbDim+1),Gbest(1:ProbDim),GearRatio);
    disp(X);
end
% g1 = 85.334407 + 0.0056858*Gbest( 2)*Gbest( 5) + 0.0006262*Gbest( 1)*Gbest( 4) - 0.0022053*Gbest( 3)*Gbest( 5);
% g2 = 80.51249 + 0.0071317*Gbest( 2)*Gbest( 5) + 0.0029955*Gbest( 1)*Gbest( 2) - 0.0021813*Gbest( 3)^2;
% g3 = 9.300961 + 0.0047026*Gbest( 3)*Gbest( 5) + 0.0012547*Gbest( 1)*Gbest( 3) + 0.0019085*Gbest( 3)*Gbest( 4);    
% 
% X = sprintf('Best Value: %d\nBest Position: %d %d %d %d %d\nConstraints:\nG1: %d\n G2: %d\n G3: %dMean: %d\nMedian %d\nStandard Deviation: %d\nWorst Value: %d\n',PSO_Curr,Gbest(ProbDim+1),Gbest(1:ProbDim),g1,g2,g3);
% disp(X);
%movie(M,1,120);
GearRatio = (Gbest(1)*Gbest(2))/(Gbest(3)*Gbest(4));
Ans(trials,:) = [Gbest GearRatio];
end

% Get Best Fit
Vals = zeros(30,1);
for o = 1:30
    Vals(o) = Ans(o,ProbDim+1);
end

% Generate Stats
Mean = mean(Vals);
StdDev = std(Vals);
Median = median(Vals);
Worst = max(Vals);
Best = min(Vals);

% Get index of best run
BesInd = 0;
for o = 1:30
    if min(Vals) == Ans(o,ProbDim+1)
        BesInd = o;
    end
end
X = sprintf('\n\nBest OverAll Value: %0.15f\nPosition: %0.15f %0.15f %0.15f %0.15f\nGear Ratio: %0.15f\nMean: %0.15f\nMedian: %0.15f\nStandard Deviation:%0.15f\nWorst Overall Value: %0.15f', Best, Ans(BesInd, 1), Ans(BesInd, 2), Ans(BesInd, 3), Ans(BesInd, 4), Ans(BesInd,ProbDim+2), Mean, Median, StdDev, Worst);
disp(X);