%LICNACHAN, LANCE OLIVER C.
%2014-64880

format long;
% rng('shuffle');

ProbDim = 4;
ConsNum = 7;
RunMax = 30;

convRuns=0;
Ans = zeros(RunMax, ProbDim+ConsNum+1+1);
timeRec = zeros(RunMax,1);
for trials = 1:RunMax
tic;
Y = sprintf('Trial: %d',trials);
disp(Y);
% Variables specific to the problem
DimMinMax = zeros(ProbDim, 2);
DimMinMax(1, :) = [.1 2];
DimMinMax(2, :) = [.1 10];
DimMinMax(3, :) = [.1 10];
DimMinMax(4, :) = [.1 2];

% Variables specific to the algorithm
AcceptThreshold = 1e-6;
PopNum = 20*ProbDim;
PSO_Curr = 1;
PSO_Max = 500;
c1 = 1.5;
c2 = 1.5;
wmax = 0.9;
wmin = 0.4;

TransPos = zeros(PopNum, ProbDim);
TransVel = zeros(PopNum, ProbDim);

% Initialization Step
[PosPop, VelPop] = PSO_InitPop(PopNum, ProbDim, DimMinMax);
FitVal = PSO_GetFitValues(PopNum, PosPop, ProbDim, DimMinMax);
[Pbest, Gbest]= PSO_GetInitPGBest(PopNum, PosPop, ProbDim, FitVal);

PrevDiff = 0;
while PSO_Curr <= PSO_Max
    % Evaluate
    FitVal = PSO_GetFitValues(PopNum, PosPop, ProbDim, DimMinMax);
    
    if(PSO_Curr == 1)
        PrevDiff = max(FitVal) - min(FitVal);
    else
        CurrDiff = max(FitVal) - min(FitVal);
        % disp(CurrDiff);
        % Check for population convergence
        if PrevDiff - CurrDiff < AcceptThreshold && CurrDiff < AcceptThreshold
            for i = 1:PopNum
                if min(FitVal) == FitVal(i)
                    minInd = i;
                end
                if max(FitVal) == FitVal(i)
                    maxInd = i;
                end
            end
            Tmax = 13600;
            SigMax = 30000;
            L = 14;
            P = 6000;
            E = 30e6;
            G = 12e6;
            R = sqrt(((PosPop(minInd, 2)^2)/4)+((PosPop(minInd, 1)+PosPop(minInd, 3))/2)^2);
            Sig = (6*P*L)/(PosPop(minInd, 4)*PosPop(minInd, 3)^2);
            Del = (6*P*L^3)/(E*PosPop(minInd, 4)*PosPop(minInd, 3)^3);
            Pc = ((4.013*E*sqrt(((PosPop(minInd, 3)^2)*(PosPop(minInd, 4)^6))/36))/L^2)*(1-((PosPop(minInd, 3)/(2*L))*sqrt(E/(4*G))));
            M = P*(L + (PosPop(minInd, 2)/2));
            J = 2*(sqrt(2)*PosPop(minInd, 1)*PosPop(minInd, 2)*(((PosPop(minInd, 2)^2)/4)+((PosPop(minInd, 1)+PosPop(minInd, 3))/2)^2));
            t1 = P/(sqrt(2)*PosPop(minInd, 1)*PosPop(minInd, 2));
            t2 = (M*R)/J;
            T = sqrt(t1^2 + 2*t1*t2*(PosPop(minInd, 2)/(2*R)) + t2^2);
            g1_b = T - Tmax;
            g2_b = Sig - SigMax;
            g3_b = PosPop(minInd, 1) - PosPop(minInd, 4);
            g4_b = 0.125 - PosPop(minInd, 1);
            g5_b = Del - 0.25;
            g6_b = P-Pc;
            g7_b = 0.10471*PosPop(minInd, 1)^2 + 0.04811*PosPop(minInd, 3)*PosPop(minInd, 4)*(14+PosPop(minInd, 2)) - 5;
            
            R = sqrt(((PosPop(maxInd, 2)^2)/4)+((PosPop(maxInd, 1)+PosPop(maxInd, 3))/2)^2);
            Sig = (6*P*L)/(PosPop(maxInd, 4)*PosPop(maxInd, 3)^2);
            Del = (6*P*L^3)/(E*PosPop(maxInd, 4)*PosPop(maxInd, 3)^3);
            Pc = ((4.013*E*sqrt(((PosPop(maxInd, 3)^2)*(PosPop(maxInd, 4)^6))/36))/L^2)*(1-((PosPop(maxInd, 3)/(2*L))*sqrt(E/(4*G))));
            M = P*(L + (PosPop(maxInd, 2)/2));
            J = 2*(sqrt(2)*PosPop(maxInd, 1)*PosPop(maxInd, 2)*(((PosPop(maxInd, 2)^2)/4)+((PosPop(maxInd, 1)+PosPop(maxInd, 3))/2)^2));
            t1 = P/(sqrt(2)*PosPop(maxInd, 1)*PosPop(maxInd, 2));
            t2 = (M*R)/J;
            T = sqrt(t1^2 + 2*t1*t2*(PosPop(maxInd, 2)/(2*R)) + t2^2);
            g1_w = T - Tmax;
            g2_w = Sig - SigMax;
            g3_w = PosPop(maxInd, 1) - PosPop(maxInd, 4);
            g4_w = 0.125 - PosPop(maxInd, 1);
            g5_w = Del - 0.25;
            g6_w = P-Pc;
            g7_w = 0.10471*PosPop(maxInd, 1)^2 + 0.04811*PosPop(maxInd, 3)*PosPop(maxInd, 4)*(14+PosPop(maxInd, 2)) - 5;
            
            mPos = mean(FitVal);
            X = sprintf('Population Converged!\nNumber of Iterations: %d\nBest Value: %0.15f\nBest Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nG5: %0.15f\nG6: %0.15f\nG7: %0.15f\nWorst Value: %0.15f\nWorst Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nG5: %0.15f\nG6: %0.15f\nG7: %0.15f\nMean: %0.15e\n',PSO_Curr,FitVal(minInd),PosPop(minInd, :),g1_b,g2_b,g3_b,g4_b,g5_b,g6_b,g7_b,FitVal(maxInd),PosPop(maxInd, :),g1_w,g2_w,g3_w,g4_w,g5_w,g6_w,g7_w,mPos);
            disp(X);
            convRuns = convRuns + 1;
            break;
        end
        PrevDiff = CurrDiff;
    end
    
    if PSO_Curr == PSO_Max
        break;
    end
    
    % Get best values
    [Pbest, Gbest] = PSO_GetPGBest(PopNum, PosPop, ProbDim, FitVal, Pbest, Gbest);
    
    % Change value according to how current iteration
    w = wmax-(wmax-wmin)*(PSO_Curr/PSO_Max);
    
    % Calculate new velocities and move
    [TransPos, TransVel] = PSO_ChangeVel(PopNum, PosPop, VelPop, ProbDim, Pbest, Gbest, w, c1, c2, DimMinMax);
    
    PosPop = TransPos;
    VelPop = TransVel;
    PSO_Curr = PSO_Curr + 1;
end

if PSO_Curr >= PSO_Max
    for i = 1:PopNum
        if min(FitVal) == FitVal(i)
            minInd = i;
        end
        if max(FitVal) == FitVal(i)
            maxInd = i;
        end
    end
    Tmax = 13600;
    SigMax = 30000;
    L = 14;
    P = 6000;
    E = 30e6;
    G = 12e6;
    R = sqrt(((PosPop(minInd, 2)^2)/4)+((PosPop(minInd, 1)+PosPop(minInd, 3))/2)^2);
    Sig = (6*P*L)/(PosPop(minInd, 4)*PosPop(minInd, 3)^2);
    Del = (6*P*L^3)/(E*PosPop(minInd, 4)*PosPop(minInd, 3)^3);
    Pc = ((4.013*E*sqrt(((PosPop(minInd, 3)^2)*(PosPop(minInd, 4)^6))/36))/L^2)*(1-((PosPop(minInd, 3)/(2*L))*sqrt(E/(4*G))));
    M = P*(L + (PosPop(minInd, 2)/2));
    J = 2*(sqrt(2)*PosPop(minInd, 1)*PosPop(minInd, 2)*(((PosPop(minInd, 2)^2)/4)+((PosPop(minInd, 1)+PosPop(minInd, 3))/2)^2));
    t1 = P/(sqrt(2)*PosPop(minInd, 1)*PosPop(minInd, 2));
    t2 = (M*R)/J;
    T = sqrt(t1^2 + 2*t1*t2*(PosPop(minInd, 2)/(2*R)) + t2^2);
    g1_b = T - Tmax;
    g2_b = Sig - SigMax;
    g3_b = PosPop(minInd, 1) - PosPop(minInd, 4);
    g4_b = 0.125 - PosPop(minInd, 1);
    g5_b = Del - 0.25;
    g6_b = P-Pc;
    g7_b = 0.10471*PosPop(minInd, 1)^2 + 0.04811*PosPop(minInd, 3)*PosPop(minInd, 4)*(14+PosPop(minInd, 2)) - 5;

    R = sqrt(((PosPop(maxInd, 2)^2)/4)+((PosPop(maxInd, 1)+PosPop(maxInd, 3))/2)^2);
    Sig = (6*P*L)/(PosPop(maxInd, 4)*PosPop(maxInd, 3)^2);
    Del = (6*P*L^3)/(E*PosPop(maxInd, 4)*PosPop(maxInd, 3)^3);
    Pc = ((4.013*E*sqrt(((PosPop(maxInd, 3)^2)*(PosPop(maxInd, 4)^6))/36))/L^2)*(1-((PosPop(maxInd, 3)/(2*L))*sqrt(E/(4*G))));
    M = P*(L + (PosPop(maxInd, 2)/2));
    J = 2*(sqrt(2)*PosPop(maxInd, 1)*PosPop(maxInd, 2)*(((PosPop(maxInd, 2)^2)/4)+((PosPop(maxInd, 1)+PosPop(maxInd, 3))/2)^2));
    t1 = P/(sqrt(2)*PosPop(maxInd, 1)*PosPop(maxInd, 2));
    t2 = (M*R)/J;
    T = sqrt(t1^2 + 2*t1*t2*(PosPop(maxInd, 2)/(2*R)) + t2^2);
    g1_w = T - Tmax;
    g2_w = Sig - SigMax;
    g3_w = PosPop(maxInd, 1) - PosPop(maxInd, 4);
    g4_w = 0.125 - PosPop(maxInd, 1);
    g5_w = Del - 0.25;
    g6_w = P-Pc;
    g7_w = 0.10471*PosPop(maxInd, 1)^2 + 0.04811*PosPop(maxInd, 3)*PosPop(maxInd, 4)*(14+PosPop(maxInd, 2)) - 5;

    mPos = mean(FitVal);
    X = sprintf('Did Not Converge!\nNumber of Iterations: %d\nBest Value: %0.15f\nBest Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nG5: %0.15f\nG6: %0.15f\nG7: %0.15f\nWorst Value: %0.15f\nWorst Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nG5: %0.15f\nG6: %0.15f\nG7: %0.15f\nMean: %0.15e\n',PSO_Curr,FitVal(minInd),PosPop(minInd, :),g1_b,g2_b,g3_b,g4_b,g5_b,g6_b,g7_b,FitVal(maxInd),PosPop(maxInd, :),g1_w,g2_w,g3_w,g4_w,g5_w,g6_w,g7_w,mPos);
    disp(X);
end
    if PSO_Curr >= PSO_Max
        Ans(trials,:) = [PosPop(minInd, :) FitVal(minInd) g1_b g2_b g3_b g4_b g5_b g6_b g7_b 0];
    else % Converged
        Ans(trials,:) = [PosPop(minInd, :) FitVal(minInd) g1_b g2_b g3_b g4_b g5_b g6_b g7_b 1];
    end
    
    timeRec(trials) = toc;
    X = sprintf('Running Time for this trial: %0.15e\n', timeRec(trials));
    disp(X);
end

% Get Best Fit
Vals = zeros(RunMax,1);
for o = 1:RunMax
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
minConvVal = 0;
for o = 1:RunMax
    if min(Vals) == Ans(o,ProbDim+1)
        BesInd = o;
    end
end

if convRuns > 0
    % Get Best Fit
    ConvVals = zeros(convRuns,1);
    i=1;
    for o = 1:RunMax
        if(Ans(o,ProbDim+ConsNum+1+1) == 1)
            ConvVals(i) = Ans(o,ProbDim+1);
            i=i+1;
        end
    end

    Best = min(ConvVals);
    for o = 1:RunMax
        if min(ConvVals) == Ans(o,ProbDim+1)
            BesInd = o;
        end
    end
end

ConvRatio = convRuns/RunMax;
totalTime = sum(timeRec);
aveTime = mean(timeRec);

disp(BesInd);
X = sprintf('\n\nBest OverAll Value: %0.15f\nPosition: %0.15f %0.15f %0.15f %0.15f\nConstraints:\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nG5: %0.15f\nG6: %0.15f\nG7: %0.15f\nMean: %0.15f\nMedian: %0.15f\nStandard Deviation:%0.15f\nWorst Best Overall Value: %0.15f\nNumber of Converged Runs: %d\nRatio of Convergence: %0.15e\nTotal Running Time for all trials: %0.15e\nAverage running time: %0.15e\n', Best, Ans(BesInd, 1), Ans(BesInd, 2), Ans(BesInd, 3), Ans(BesInd, 4), Ans(BesInd, 6), Ans(BesInd, 7), Ans(BesInd, 8), Ans(BesInd, 9),  Ans(BesInd, 10), Ans(BesInd, 11), Ans(BesInd, 12), Mean, Median, StdDev, Worst, convRuns, ConvRatio, totalTime, aveTime);
disp(X);
