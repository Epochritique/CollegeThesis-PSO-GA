%LICNACHAN, LANCE OLIVER C.
%2014-64880

ProbDim = 5;

Ans = zeros(30, ProbDim+4);
for CD = 1:30
X = sprintf('Iteration: %d\n',CD);
disp(X);
format long;

ctr=0;

% Variables specific to the problem
DimMinMax = zeros(ProbDim, 2);
DimMinMax(1, :) = [78 102];
DimMinMax(2, :) = [33 45];
DimMinMax(3, :) = [27 45];
DimMinMax(4, :) = [27 45];
DimMinMax(5, :) = [27 45];

InitVeloMax = -10;
InitVeloMin = 10;

% Variables specific to the algorithm
AcceptThreshold = 1e-6;
Aim = -30373.9490;

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
       TempPos(1, i) = (DimMinMax(i, 2)- DimMinMax(i, 1))*rand(1,1,'double') + DimMinMax(i, 1);
    end
    PosPop(RowNum, :) = TempPos; % Add to positions
    TempVel = zeros(1, ProbDim);
    for i = 1:ProbDim
       TempVel(1, i) = (InitVeloMax-InitVeloMin)*rand(1,1,'double') + InitVeloMin;
    end
    VelPop(RowNum, :) = TempVel; % Add to velocities
end

% Evaluate each particle's value
for RowNum = 1:PopNum
    FitVal(RowNum, 1) = 5.3578547*PosPop(RowNum, 3)^2 + 0.8356891*PosPop(RowNum, 1)*PosPop(RowNum, 5) + 37.293239*PosPop(RowNum, 1) - 40792.141;
end
FitWorst = max(FitVal);
for RowNum = 1:PopNum
    probs = 0;
    for i = 1:ProbDim
       if ~ (PosPop(RowNum, i)>= DimMinMax(i, 1) && PosPop(RowNum, i)<= DimMinMax(i, 2))
           probs=probs+1;% compute the number of errors
       end
    end

    % Check if it satisfies constraint equations
    g1 = 85.334407 + 0.0056858*PosPop(RowNum, 2)*PosPop(RowNum, 5) + 0.0006262*PosPop(RowNum, 1)*PosPop(RowNum, 4) - 0.0022053*PosPop(RowNum, 3)*PosPop(RowNum, 5);
    g2 = 80.51249 + 0.0071317*PosPop(RowNum, 2)*PosPop(RowNum, 5) + 0.0029955*PosPop(RowNum, 1)*PosPop(RowNum, 2) - 0.0021813*PosPop(RowNum, 3)^2;
    g3 = 9.300961 + 0.0047026*PosPop(RowNum, 3)*PosPop(RowNum, 5) + 0.0012547*PosPop(RowNum, 1)*PosPop(RowNum, 3) + 0.0019085*PosPop(RowNum, 3)*PosPop(RowNum, 4);
    if ~(g1>=0 && g1<=92)
            probs=probs+1;% compute the number of errors
    end
    if ~(g2>=90 && g2<=110)
            probs=probs+1;% compute the number of errors
    end
    if ~(g3>=20 && g3<=25)
            probs=probs+1;% compute the number of errors
    end
    if(probs~=0)
        FitVal(RowNum, 1) = FitWorst + g1 + g2 + g3;
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
        FitVal(RowNum, 1) = 5.3578547*PosPop(RowNum, 3)^2 + 0.8356891*PosPop(RowNum, 1)*PosPop(RowNum, 5) + 37.293239*PosPop(RowNum, 1) - 40792.141;
    end
    FitWorst = max(FitVal);
    for RowNum = 1:PopNum
        probs = 0;
        for i = 1:ProbDim
           if ~ (PosPop(RowNum, i)>= DimMinMax(i, 1) && PosPop(RowNum, i)<= DimMinMax(i, 2))
               probs=probs+1;% compute the number of errors
           end
        end

        % Check if it satisfies constraint equations
        g1 = 85.334407 + 0.0056858*PosPop(RowNum, 2)*PosPop(RowNum, 5) + 0.0006262*PosPop(RowNum, 1)*PosPop(RowNum, 4) - 0.0022053*PosPop(RowNum, 3)*PosPop(RowNum, 5);
        g2 = 80.51249 + 0.0071317*PosPop(RowNum, 2)*PosPop(RowNum, 5) + 0.0029955*PosPop(RowNum, 1)*PosPop(RowNum, 2) - 0.0021813*PosPop(RowNum, 3)^2;
        g3 = 9.300961 + 0.0047026*PosPop(RowNum, 3)*PosPop(RowNum, 5) + 0.0012547*PosPop(RowNum, 1)*PosPop(RowNum, 3) + 0.0019085*PosPop(RowNum, 3)*PosPop(RowNum, 4);
        if ~(g1>=0 && g1<=92)
                probs=probs+1;% compute the number of errors
        end
        if ~(g2>=90 && g2<=110)
                probs=probs+1;% compute the number of errors
        end
        if ~(g3>=20 && g3<=25)
                probs=probs+1;% compute the number of errors
        end
        if(probs~=0)
            FitVal(RowNum, 1) = FitWorst + g1 + g2 + g3;
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
        TempVel = VelPop(RowNum, 1:ProbDim)*w + (c1*rand())*(Pbest(RowNum, 1:ProbDim) - PosPop(RowNum, 1:ProbDim)) + (c2*rand())*(Gbest(1:ProbDim) - PosPop(RowNum, 1:ProbDim));
        TempPos = PosPop(RowNum, 1) + TempVel;
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
                   TempPos(1, i) = (DimMinMax(i, 2) - DimMinMax(i, 1))*rand(1,1,'double') + DimMinMax(i, 1);
                end
                GA_Chroms(RowNum, :) = TempPos; % Add to positions
            end
            
            % Set the first chromosome as the selected chromosome
            GA_Chroms(1, :) = Sel_Indiv;

            % Compute Fitness Values
            GA_FitVal = zeros(GA_PS, 1);
            for RowNum = 1:GA_PS
                GA_FitVal(RowNum, 1) = 5.3578547*GA_Chroms(RowNum, 3)^2 + 0.8356891*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 5) + 37.293239*GA_Chroms(RowNum, 1) - 40792.141;
            end
            GA_FitWorst = max(GA_FitVal);
            for RowNum = 1:GA_PS
                probs = 0;
                for i = 1:ProbDim
                   if ~ (GA_Chroms(RowNum, i)>= DimMinMax(i, 1) && GA_Chroms(RowNum, i)<= DimMinMax(i, 2))
                       probs=probs+1;% compute the number of errors
                   end
                end

                % Check if it satisfies constraint equations
                GA_g1 = 85.334407 + 0.0056858*GA_Chroms(RowNum, 2)*GA_Chroms(RowNum, 5) + 0.0006262*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 4) - 0.0022053*GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 5);
                GA_g2 = 80.51249 + 0.0071317*GA_Chroms(RowNum, 2)*GA_Chroms(RowNum, 5) + 0.0029955*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 2) - 0.0021813*GA_Chroms(RowNum, 3)^2;
                GA_g3 = 9.300961 + 0.0047026*GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 5) + 0.0012547*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 3) + 0.0019085*GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 4);
                if ~(GA_g1>=0 && GA_g1<=92)
                        probs=probs+1;% compute the number of errors
                end
                if ~(GA_g2>=90 && GA_g2<=110)
                        probs=probs+1;% compute the number of errors
                end
                if ~(GA_g3>=20 && GA_g3<=25)
                        probs=probs+1;% cofmpute the number of errors
                end
                if(probs~=0)
                    GA_FitVal(RowNum, 1) = GA_FitWorst + GA_g1 + GA_g2 + GA_g3;
                end
            end

            % Start GA process
            GA_Curr=0;
            while GA_Curr <= GA_MaxItr
                % Build Roulette Wheel for GA Selection 
                GA_FitRatio = (1/((GA_FitVal-Aim)+1)) / sum(1/((GA_FitVal-Aim)+1)); 
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

                    RanSelVal = sum(1/(GA_FitVal-Aim+1))*rand(1,'double');
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

                    RanSelVal = sum(1/(GA_FitVal-Aim+1))*rand(1,'double');
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
                    Sib1_FitVal = 5.3578547*Sibling1(3)^2 + 0.8356891*Sibling1(1)*Sibling1(5) + 37.293239*Sibling1(1) - 40792.141;
                    Sib2_FitVal = 5.3578547*Sibling2(3)^2 + 0.8356891*Sibling2(1)*Sibling2(5) + 37.293239*Sibling2(1) - 40792.141;
                    if Sib1_FitVal >= Sib2_FitVal
                        SibRep = Sibling1;
                    else
                        SibRep = Sibling2;
                    end

                    % Mutation
                    MutVal = rand();
                    if MutVal <= GA_mut
                        MutPos = randi(ProbDim);
                        SibRep(1,MutPos) = (DimMinMax(MutPos, 2)- DimMinMax(MutPos, 1))*rand(1,'double') + DimMinMax(MutPos, 1);
                    end
                    GA_TransPop(Replacements, :) = SibRep(1, :);
                    Replacements = Replacements + 1;
                end
                GA_Chroms = GA_TransPop;
                
                % Compute Fitness Values
                GA_FitVal = zeros(GA_PS, 1);
                for RowNum = 1:GA_PS
                    GA_FitVal(RowNum, 1) = 5.3578547*GA_Chroms(RowNum, 3)^2 + 0.8356891*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 5) + 37.293239*GA_Chroms(RowNum, 1) - 40792.141;
                end
                GA_FitWorst = max(GA_FitVal);

                for RowNum = 1:GA_PS
                    probs = 0;
                    for i = 1:ProbDim
                       if ~ (GA_Chroms(RowNum, i)>= DimMinMax(i, 1) && GA_Chroms(RowNum, i)<= DimMinMax(i, 2))
                           probs=probs+1;% compute the number of errors
                       end
                    end

                    % Check if it satisfies constraint equations
                    GA_g1 = 85.334407 + 0.0056858*GA_Chroms(RowNum, 2)*GA_Chroms(RowNum, 5) + 0.0006262*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 4) - 0.0022053*GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 5);
                    GA_g2 = 80.51249 + 0.0071317*GA_Chroms(RowNum, 2)*GA_Chroms(RowNum, 5) + 0.0029955*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 2) - 0.0021813*GA_Chroms(RowNum, 3)^2;
                    GA_g3 = 9.300961 + 0.0047026*GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 5) + 0.0012547*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 3) + 0.0019085*GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 4);
                    if ~(GA_g1>=0 && GA_g1<=92)
                            probs=probs+1;% compute the number of errors
                    end
                    if ~(GA_g2>=90 && GA_g2<=110)
                            probs=probs+1;% compute the number of errors
                    end
                    if ~(GA_g3>=20 && GA_g3<=25)
                            probs=probs+1;% compute the number of errors
                    end
                    if(probs~=0)
                        GA_FitVal(RowNum, 1) = GA_FitWorst + GA_g1 + GA_g2 + GA_g3;
                    end
                end

                Best_GA = zeros(1, ProbDim);
                for RowNum = 1:GA_PS
                    if min(GA_FitVal) == GA_FitVal(RowNum)
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
    
    if (Gbest(ProbDim+1) <= Aim)
        ctr = ctr+1;
        if ctr >= 5
            X = sprintf('Converged!\nNumber of Iterations: %d\nBest Value: %d\nBest Location: %d %d %d %d %d\n',PSO_Curr,Gbest(ProbDim+1),Gbest(1), Gbest(2), Gbest(3), Gbest(4), Gbest(5));
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
    X = sprintf('Did Not Converge!\nNumber of Iterations: %d\nBest Value: %d\nBest Location: %d %d %d %d %d\n',PSO_Curr,Gbest(ProbDim+1),Gbest(1), Gbest(2), Gbest(3), Gbest(4), Gbest(5));
    disp(X);
end

g1 = 85.334407 + 0.0056858*Gbest(2)*Gbest(5) + 0.0006262*Gbest(1)*Gbest(4) - 0.0022053*Gbest(3)*Gbest(5);
g2 = 80.51249 + 0.0071317*Gbest(2)*Gbest(5) + 0.0029955*Gbest(1)*Gbest(2) - 0.0021813*Gbest(3)^2;
g3 = 9.300961 + 0.0047026*Gbest(3)*Gbest(5) + 0.0012547*Gbest(1)*Gbest(3) + 0.0019085*Gbest(3)*Gbest(4);
Ans(CD,:) = [Gbest g1 g2 g3];

end

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

X = sprintf('\n\nBest OverAll Value: %d\nPosition: %d %d %d %d %d\nConstraints:\nG1: %d\nG2: %d\nG3: %d\nMean: %d\nMedian: %d\nStandard Deviation:%d\nWorst Overall Valu: %d', Best, Ans(BesInd, 1), Ans(BesInd, 2), Ans(BesInd, 3), Ans(BesInd, 4), Ans(BesInd, 5), Ans(BesInd, 7), Ans(BesInd, 8), Ans(BesInd, 9), Mean, Median, StdDev, Worst);
disp(X);