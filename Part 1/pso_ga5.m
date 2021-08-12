%LICNACHAN, LANCE OLIVER C.
%2014-64880

ProbDim = 4;
ConsProbs = 0;

Ans = zeros(30, ProbDim+ConsProbs+3);
for CD = 1:30
X = sprintf('Iteration: %d\n',CD);
disp(X);
format long;

ctr=0;

% Variables specific to the problem
DimMinMax = zeros(ProbDim, 2);
DimMinMax(1, :) = [12 60];
DimMinMax(2, :) = [12 60];
DimMinMax(3, :) = [12 60];
DimMinMax(4, :) = [12 60];

InitVeloMax = 10;
InitVeloMin = -10;

% Variables specific to the algorithm
AcceptThreshold = 1e-6;
Aim = 2.70085e-12;
PopNum = 20;

PSO_Curr=1;
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
GA_MaxItr = GA_MinItr + ((PSO_Curr/PSO_Max)^GA_B)*(GA_MaxItr-GA_MinItr);
GA_Num = GA_NumMax - ((PSO_Curr/PSO_Max)^GA_y)*(GA_NumMax-GA_NumMin);
GA_PS = ceil(GA_MinPS + ((PSO_Curr/PSO_Max)^GA_y)*(GA_MaxPS-GA_MinPS));

PosPop = zeros(PopNum, ProbDim);
VelPop = zeros(PopNum, ProbDim);
TransPos = zeros(PopNum, ProbDim);
TransVel = zeros(PopNum, ProbDim);
RowNum = 1;
FitVal = zeros(PopNum, 1);
Pbest = zeros(PopNum, ProbDim+1);
Gbest = zeros(1, ProbDim+1);

% Get initial population
for RowNum = 1:PopNum
    TempPos = zeros(1, ProbDim);
    for i = 1:ProbDim
       TempPos(1, i) = randi([DimMinMax(i, 1), DimMinMax(i, 2)], 1);
    end
    PosPop(RowNum, :) = TempPos; % Add to positions
    TempVel = zeros(1, ProbDim);
    for i = 1:ProbDim
    	TempVel(1, i) = randi([InitVeloMin, InitVeloMax], 1);
    end
    VelPop(RowNum, :) = TempVel; % Add to velocities
end

% Evaluate each particle's value
for RowNum = 1:PopNum
    FitVal(RowNum, 1) = ((1/6.931) - ((PosPop(RowNum, 1)*PosPop(RowNum, 2))/(PosPop(RowNum, 3)*PosPop(RowNum, 4))))^2;
end
FitWorst = max(FitVal);
for RowNum = 1:PopNum
    probs = 0;
    for i = 1:ProbDim
       if ~ (PosPop(RowNum, i)>= DimMinMax(i, 1) && PosPop(RowNum, i)<= DimMinMax(i, 2))
           probs=probs+1;% compute the number of errors
       end
    end
    if(probs~=0)
        FitVal(RowNum, 1) = abs(FitWorst);
    end
end

% Get Initial Pbest and Gbest
for RowNum = 1:PopNum
    Pbest(RowNum, :) = [PosPop(RowNum,:) FitVal(RowNum, 1)];
end

Gbest(ProbDim+1) = min(FitVal);
for RowNum = 1:PopNum
    if Pbest(RowNum, ProbDim + 1) == Gbest(ProbDim+1)
        Gbest = Pbest(RowNum, :);
        break;
    end
end

while PSO_Curr <= PSO_Max
    % Evaluate each particle's value
    for RowNum = 1:PopNum
        FitVal(RowNum, 1) = ((1/6.931) - ((PosPop(RowNum, 1)*PosPop(RowNum, 2))/(PosPop(RowNum, 3)*PosPop(RowNum, 4))))^2;
    end
    FitWorst = max(FitVal);
    for RowNum = 1:PopNum
        probs = 0;
        for i = 1:ProbDim
           if ~ (PosPop(RowNum, i)>= DimMinMax(i, 1) && PosPop(RowNum, i)<= DimMinMax(i, 2))
               probs=probs+1;% compute the number of errors
           end
        end
        if(probs~=0)
            FitVal(RowNum, 1) = abs(FitWorst);
        end
    end


    % Get Pbest and Gbest
    for RowNum = 1:PopNum
        if Pbest(RowNum, ProbDim+1) > FitVal(RowNum,1)
            Pbest(RowNum, :) = [PosPop(RowNum,:) FitVal(RowNum, 1)];
        end
    end
    
    if Gbest(ProbDim+1) > min(FitVal)
        Gbest(ProbDim+1) = min(FitVal);
        for RowNum = 1:PopNum
            if Pbest(RowNum, ProbDim + 1) == Gbest(ProbDim+1)
                Gbest = Pbest(RowNum, :);
                break;
            end
        end
    end
    
    % Change Velocity According to Best Positions
    w = wmax-(wmax-wmin)*(PSO_Curr/PSO_Max);

    for RowNum = 1:PopNum
        TempVel = ceil(VelPop(RowNum, 1:ProbDim)*w + (c1*rand())*(Pbest(RowNum, 1:ProbDim) - PosPop(RowNum, 1:ProbDim)) + (c2*rand())*(Gbest(1:ProbDim) - PosPop(RowNum, 1:ProbDim)));
        TempPos = ceil(PosPop(RowNum, 1) + TempVel);
        TransPos(RowNum, :) = TempPos;
        TransVel(RowNum, :) = TempVel;
    end
    
        PSO_Arranged = sort(FitVal);
        
        % GA Portion
        GA_Num_Curr = 1;
        while GA_Num_Curr <= GA_Num
            % Get one from best individuals
            for RowNum = 1:PopNum
                if FitVal(RowNum) == PSO_Arranged(GA_Num_Curr);
                   Sel_Indiv = PosPop(RowNum, :); 
                end
            end

            GA_Chroms = zeros(GA_PS, ProbDim);
            % Randomly Generate GA_PS chromosomes
            for RowNum = 2:GA_PS
                TempPos = zeros(1, ProbDim);
                for i = 1:ProbDim
                   TempPos(1, i) = randi([DimMinMax(i, 1), DimMinMax(i, 2)], 1);
                end
                GA_Chroms(RowNum, :) = TempPos; % Add to positions
            end
            
            % Set the first chromosome as the selected chromosome
            GA_Chroms(1, :) = Sel_Indiv;

            % Compute Fitness Values
            GA_FitVal = zeros(GA_PS, 1);
            for RowNum = 1:GA_PS
                GA_FitVal(RowNum, 1) = ((1/6.931) - ((GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 2))/(GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 4))))^2;
            end
            GA_FitWorst = max(GA_FitVal);
            for RowNum = 1:GA_PS
                probs = 0;
                for i = 1:ProbDim
                   if ~ (GA_Chroms(RowNum, i)>= DimMinMax(i, 1) && GA_Chroms(RowNum, i)<= DimMinMax(i, 2))
                       probs=probs+1;% compute the number of errors
                   end
                end
                if(probs~=0)
                    GA_FitVal(RowNum, 1) = abs(GA_FitWorst);
                end
            end

            % Start GA process
            GA_Curr=0;
            while GA_Curr <= GA_MaxItr
                % Build Roulette Wheel for GA Selection 
                GA_FitRatio = abs(log(GA_FitVal)) / sum(abs(log(GA_FitVal))); 
                GA_RouWheel = zeros(GA_PS, 1);
                for RowNum = 1:GA_PS
                    for i = 1:RowNum
                        GA_RouWheel(RowNum) = GA_RouWheel(RowNum) + GA_FitRatio(i);
                    end
                end

                % Transition Population with Elitist Mode
                GA_TransPop = zeros(GA_PS, ProbDim);
                Arranged = sort(GA_FitVal);
                for RowNum = 1:floor(GA_PS/2)
                    for i = 1:GA_PS
                        if GA_FitVal(i) == Arranged(RowNum)
                            GA_TransPop(RowNum, :) = GA_Chroms(i, :);
                        end
                    end
                end

                Replacements = floor(GA_PS/2) + 1;
                while Replacements <= GA_PS
                    % Select 2 for Cross_Over
                    Parent1 = zeros(1, ProbDim);
                    Parent2 = zeros(1, ProbDim);

                    RanSelVal = sum(abs(log(GA_FitVal)))*rand(1,'double');
                    for RowNum = 1:GA_PS
                        if GA_RouWheel(RowNum) > RanSelVal
                            if RowNum == 1
                                Parent1(1, :) = GA_Chroms(RowNum, :);
                            else 
                                Parent1(1, :) = GA_Chroms(RowNum-1, :);
                            end
                            break;
                        end
                    end

                    RanSelVal = sum(abs(log(GA_FitVal)))*rand(1,'double');
                    for RowNum = 1:GA_PS
                        if GA_RouWheel(RowNum) > RanSelVal
                            if RowNum == 1
                                Parent2(1, :) = GA_Chroms(RowNum, :);
                            else 
                                Parent2(1, :) = GA_Chroms(RowNum-1, :);
                            end 
                            break;
                        end
                    end
                    
                    % Cross_Over
                    CrossVal = rand();
                    if CrossVal <= GA_cross
                        CrossPos = randi(ProbDim-1);
                        Sibling1 = [Parent1(1:CrossPos) Parent2(CrossPos+1:ProbDim)];
                        Sibling2 = [Parent2(1:CrossPos) Parent1(CrossPos+1:ProbDim)];
                    else
                        Sibling1 = Parent1;
                        Sibling2 = Parent2;
                    end

                    % Compute Fitness Values
                    Sib1_FitVal = ((1/6.931) - ((Sibling1(1)*Sibling1(2))/(Sibling1(3)*Sibling1(4))))^2;
                    Sib2_FitVal = ((1/6.931) - ((Sibling2(1)*Sibling2(2))/(Sibling2(3)*Sibling2(4))))^2;
                    if Sib1_FitVal >= Sib2_FitVal
                        SibRep = Sibling1;
                    else
                        SibRep = Sibling2;
                    end

                    % Mutation
                    MutVal = rand();
                    if MutVal <= GA_mut
                        MutPos = randi(ProbDim);
                        SibRep(1,MutPos) = randi([DimMinMax(MutPos, 1), DimMinMax(MutPos, 2)],1);
                    end
                    GA_TransPop(Replacements, :) = SibRep(1, :);
                    Replacements = Replacements + 1;
                end
                GA_Chroms = GA_TransPop;
                
                % Compute Fitness Values
                GA_FitVal = zeros(GA_PS, 1);
                for RowNum = 1:GA_PS
                    GA_FitVal(RowNum, 1) = ((1/6.931) - ((GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 2))/(GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 4))))^2;
                end
                GA_FitWorst = max(GA_FitVal);
                for RowNum = 1:GA_PS
                    probs = 0;
                    for i = 1:ProbDim
                       if ~ (GA_Chroms(RowNum, i)>= DimMinMax(i, 1) && GA_Chroms(RowNum, i)<= DimMinMax(i, 2))
                           probs=probs+1;% compute the number of errors
                       end
                    end
                    if(probs~=0)
                        GA_FitVal(RowNum, 1) = abs(GA_FitWorst);
                    end
                end

                Best_GA = zeros(1, ProbDim);
                for RowNum = 1:GA_PS
                    if min(GA_FitVal) == FitVal(RowNum)
                        Best_GA(1, :) = GA_Chroms(RowNum, :);
                    end
                end
                GA_Curr = GA_Curr + 1;
            end
            
            % Replace worst individual
            for RowNum = 1:PopNum
                if FitVal(RowNum) == PSO_Arranged((ceil(GA_Num+1))-GA_Num_Curr);
                   PosPop(RowNum,:) = Best_GA(1,:);
                end
            end
            GA_Num_Curr = GA_Num_Curr + 1;
        end
        % Update GA_Vars
        GA_MaxItr = GA_MinItr + ((PSO_Curr/PSO_Max)^GA_B)*(GA_MaxItr-GA_MinItr);
        GA_Num = ceil(GA_NumMax - ((PSO_Curr/PSO_Max)^GA_y)*(GA_NumMax-GA_NumMin));
        GA_PS = ceil(GA_MinPS + ((PSO_Curr/PSO_Max)^GA_y)*(GA_MaxPS-GA_MinPS));

    GearRatio = (Gbest(1)*Gbest(2))/(Gbest(3)*Gbest(4));
    if Gbest(ProbDim+1) <= Aim
        ctr = ctr+1;
        if ctr >= 5
            X = sprintf('Converged!\nNumber of Iterations: %d\nBest Value: %d\nBest Location: %d %d %d %d\nGear Ratio: %d',PSO_Curr,Gbest(ProbDim+1),Gbest(1), Gbest(2), Gbest(3), Gbest(4), GearRatio);
            disp(X);
            break;
        end
    else
        ctr = 0;
    end
    
    PosPop = TransPos;
    VelPop = TransVel;
    PSO_Curr=PSO_Curr+1;
end

if PSO_Curr > PSO_Max
    X = sprintf('Did Not Converge!\nNumber of Iterations: %d\nBest Value: %d\nBest Location: %d %d %d %d\nGear Ratio: %d',PSO_Curr,Gbest(ProbDim+1),Gbest(1), Gbest(2), Gbest(3), Gbest(4), GearRatio);
    disp(X);
end

if PSO_Curr > PSO_Max
    conv = 0;
else
    conv = 1;
end
Ans(CD,:) = [Gbest GearRatio conv];
end

convs = 0;
for o = 1:30
    if Ans(o, ProbDim+ConsProbs+3) == 1
        convs = convs + 1;
    end
end

if convs == 0
    Vals = zeros(30,1);
    for o = 1:30
        Vals(o) = Ans(o,ProbDim+1);
    end

    % Mean
    Mean = mean(Vals);
    StdDev = std(Vals);
    Median = median(Vals);
    Worst = max(Vals);
    Best = min(Vals);

    BesInd = 0;
    for o = 1:30
        if min(Vals) == Ans(o,ProbDim+1)
            BesInd = o;
        end
    end

    X = sprintf('\n\nBest OverAll Value: %d\nPosition: %d %d %d %d\nMean: %d\nMedian: %d\nStandard Deviation:%d\nWorst Overall Valu: %d', Best, Ans(BesInd, 1), Ans(BesInd, 2), Ans(BesInd, 3), Ans(BesInd, 4), Mean, Median, StdDev, Worst);
    disp(X);
else
    i = 1;
    Vals = zeros(convs,1);
    for o = 1:30
        if Ans(o, ProbDim+ConsProbs+3) == 1
            Vals(i) = Ans(o,ProbDim+1);
            i = i+1;
        end
    end

    % Mean
    Mean = mean(Vals);
    StdDev = std(Vals);
    Median = median(Vals);
    Worst = max(Vals);
    Best = min(Vals);

    BesInd = 0;
    for o = 1:30
        if min(Vals) == Ans(o,ProbDim+1)
            BesInd = o;
        end
    end

    X = sprintf('\n\nBest OverAll Value: %d\nPosition: %d %d %d %d\nGear Ratio: %d\nMean: %d\nMedian: %d\nStandard Deviation:%d\nWorst Overall Valu: %d', Best, Ans(BesInd, 1), Ans(BesInd, 2), Ans(BesInd, 3), Ans(BesInd, 4), Ans(BesInd,ProbDim+2), Mean, Median, StdDev, Worst);
    disp(X);
end