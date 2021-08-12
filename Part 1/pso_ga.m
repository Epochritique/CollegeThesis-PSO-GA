%LICNACHAN, LANCE OLIVER C.
%2014-64880

format long;

% Variables specific to the problem
ProbDim = 5;
DimMinMax = zeros(ProbDim, 2);
DimMinMax(1, :) = [78 102];
DimMinMax(2, :) = [33 45];
DimMinMax(3, :) = [27 45];
DimMinMax(4, :) = [27 45];
DimMinMax(5, :) = [27 45];

InitVeloMax = -15;
InitVeloMin = 15;

% Variables specific to the algorithm
AcceptThreshold = 1e-6;
PopNum = 20*ProbDim;

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
GA_MaxItr = 0;
GA_MaxItr = GA_MinItr + ((PSO_Curr/PSO_Max)^GA_B)*(GA_MaxItr-GA_MinItr);
GA_Num = GA_NumMax - ((PSO_Curr/PSO_Max)^GA_y)*(GA_NumMax-GA_NumMin);
GA_PS = GA_MinPS + ((PSO_Curr/PSO_Max)^GA_y)*(GA_MaxPS-GA_MinPS);

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
       TempPos(1, i) = (DimMinMax(i, 2) - DimMinMax(i, 1))*rand(1,'double') + DimMinMax(i, 1);
    end
    PosPop(RowNum, :) = TempPos; % Add to positions
    TempVel = zeros(1, ProbDim);
    for i = 1:ProbDim
       TempVel(1, i) = (InitVeloMax - InitVeloMin)*rand(1,'double') + InitVeloMin;
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
       if ~ (PosPop(RowNum, i) >= DimMinMax(i, 1) && PosPop(RowNum, i) <= DimMinMax(i, 2))
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
        Gbest = [Pbest(RowNum, :) FitVal(RowNum, 1)];
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
                Gbest = [Pbest(RowNum, :) FitVal(RowNum, 1)];
                break;
            end
        end
    end
    
    % Change Velocity According to Best Positions
    w = wmax-(wmax-wmin)*(PSO_Curr/PSO_Max);

    for RowNum = 1:PopNum
        TempVel = VelPop(RowNum, 1:ProbDim)*w + c1*rand()*(Pbest(RowNum, 1:ProbDim) - PosPop(RowNum, 1:ProbDim)) + c2*rand()*(Gbest(1, 1:ProbDim) - PosPop(RowNum, 1:ProbDim));
        TempPos = PosPop(RowNum, 1) + TempVel;
        TransPos(RowNum, :) = TempPos;
        TransVel(RowNum, :) = TempVel;
    end

    
        % GA Portion
        GA_Num_Curr = 0;
        while GA_Num_Curr <= GA_Num
            GA_Num_Curr = GA_Num_Curr + 1;
            % Select a Position from PSO according to its Fitness Value
            FitRatio = abs(FitVal)/sum(abs(FitVal)); 
            
            % Build Wheel according to Fitness Value
            RouWheel = zeros(PopNum, 1);
            for RowNum = 1:PopNum
                for i = 1:RowNum
                    RouWheel(RowNum) = RouWheel(RowNum) + FitRatio(i);
                end
            end

            Sel_Indiv = zeros(1, ProbDim);
            ReplaceNumber = 1;
            % Select an individual
            RanSelVal = sum(abs(FitVal))*rand(1,1,'double');
            for RowNum = 1:PopNum
                if RouWheel(RowNum) > RanSelVal
                    if RowNum == 1
                        Sel_Indiv(1, :) = PosPop(RowNum, :);
                        ReplaceNumber = RowNum;
                    else
                        Sel_Indiv(1, :) = PosPop(RowNum-1, :);
                        ReplaceNumber = RowNum-1;
                    end
                    break;
                end
            end

            GA_Chroms = zeros(GA_PS, ProbDim);
            % Randomly Generate GA_PS chromosomes
            for RowNum = 1:GA_PS
                TempPos = zeros(1, ProbDim);
                for i = 1:ProbDim
                   TempPos(1, i) = (DimMinMax(i, 2)- DimMinMax(i, 1))*rand(1,1,'double') + DimMinMax(i, 1);
                end
                GA_Chroms(RowNum, :) = TempPos; % Add to positions
            end
            % Set the first chromosome as the selected chromosome
            GA_Chrom(1, :) = Sel_Indiv;

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
                GA_FitRatio = abs(GA_FitVal)/sum(abs(GA_FitVal)); 
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

                Replacements = floor(GA_PS/2);
                while Replacements <= GA_PS
                    % Select 2 for Cross_Over
                    Parent1 = zeros(1, ProbDim);
                    Parent2 = zeros(1, ProbDim);

                    RanSelVal = sum(abs(FitVal))*rand(1,1,'double');
                    for RowNum = 1:GA_PS
                        if GA_RouWheel(RowNum) > RanSelVal
                            if RowNum == 1
                                Parent1(1, :) = GA_Chrom(RowNum, :);
                            else 
                                Parent1(1, :) = GA_Chrom(RowNum-1, :);
                            end
                            break;
                        end
                    end

                    RanSelVal = sum(abs(FitVal))*rand(1,1,'double');
                    for RowNum = 1:GA_PS
                        if GA_RouWheel(RowNum) > RanSelVal
                            if RowNum == 1
                                Parent2(1, :) = GA_Chrom(RowNum, :);
                            else 
                                Parent2(1, :) = GA_Chrom(RowNum-1, :);
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

                    if randi(2) == 1
                        SibRep = Sibling1;
                    else
                        SibRep = Sibling2;
                    end

                    % Mutation
                    MutVal = rand();
                    if MutVal <= GA_mut
                        MutPos = randi(ProbDim);
                        SibRep(1,MutPos) = (DimMinMax(MutPos, 2)- DimMinMax(MutPos, 1))*rand(1,1,'double') + DimMinMax(MutPos, 1);
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
                    if min(GA_FitVal) == FitVal(RowNum)
                        Best_GA(1, :) = GA_Chroms(RowNum, :);
                    end
                end
                GA_Curr = GA_Curr + 1;
            end

            PosPop(ReplaceNumber,:) = Best_GA(1,:);
        end
        % Update GA_Vars
        GA_MaxItr = GA_MinItr + ((PSO_Curr/PSO_Max)^GA_B)*(GA_MaxItr-GA_MinItr);
        GA_Num = GA_NumMax - ((PSO_Curr/PSO_Max)^GA_y)*(GA_NumMax-GA_NumMin);
        GA_PS = ceil(GA_MinPS + ((PSO_Curr/PSO_Max)^GA_y)*(GA_MaxPS-GA_MinPS));

    if (Gbest(ProbDim+1) < -30600)     
        X = sprintf('Converged!\nNumber of Iterations: %d\nBest Value: %d\nBest Location: %d %d %d %d %d\n',PSO_Curr,Gbest(ProbDim+1),Gbest(1), Gbest(2), Gbest(3), Gbest(4), Gbest(5));
        disp(X);
        break;
    end
    PosPop = TransPos;
    VelPop = TransVel;
    PSO_Curr=PSO_Curr+1;
end

if PSO_Curr > PSO_Max
    X = sprintf('Did Not Converge!\nNumber of Iterations: %d\nBest Value: %d\nBest Location: %d %d %d %d %d',PSO_Curr,Gbest(ProbDim+1),Gbest(1), Gbest(2), Gbest(3), Gbest(4), Gbest(5));
    disp(X);
end
