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
GA_cross = 0.85;
GA_mut = 0.02;
GA_Curr = 1;
GA_MaxItr = 500;
GA_PS = 20*ProbDim; % Population Size

GA_Chroms = GA_InitPop(GA_PS, ProbDim, DimMinMax);

PrevDiff = 0;
while GA_Curr <= GA_MaxItr
    GA_FitVal = GA_GetFitValues(GA_PS, GA_Chroms, ProbDim, DimMinMax);
    
    if(GA_Curr == 1)
        PrevDiff = max(GA_FitVal) - min(GA_FitVal);
    else
        CurrDiff = max(GA_FitVal) - min(GA_FitVal);
        % Check for population convergence
        if PrevDiff - CurrDiff < AcceptThreshold && CurrDiff < AcceptThreshold
            for i = 1:GA_PS
                if min(GA_FitVal) == GA_FitVal(i)
                    minInd = i;
                end
                if max(GA_FitVal) == GA_FitVal(i)
                    maxInd = i;
                end
            end
            Tmax = 13600;
            SigMax = 30000;
            L = 14;
            P = 6000;
            E = 30e6;
            G = 12e6;
            R = sqrt(((GA_Chroms(minInd, 2)^2)/4)+((GA_Chroms(minInd, 1)+GA_Chroms(minInd, 3))/2)^2);
            Sig = (6*P*L)/(GA_Chroms(minInd, 4)*GA_Chroms(minInd, 3)^2);
            Del = (6*P*L^3)/(E*GA_Chroms(minInd, 4)*GA_Chroms(minInd, 3)^3);
            Pc = ((4.013*E*sqrt(((GA_Chroms(minInd, 3)^2)*(GA_Chroms(minInd, 4)^6))/36))/L^2)*(1-((GA_Chroms(minInd, 3)/(2*L))*sqrt(E/(4*G))));
            M = P*(L + (GA_Chroms(minInd, 2)/2));
            J = 2*(sqrt(2)*GA_Chroms(minInd, 1)*GA_Chroms(minInd, 2)*(((GA_Chroms(minInd, 2)^2)/4)+((GA_Chroms(minInd, 1)+GA_Chroms(minInd, 3))/2)^2));
            t1 = P/(sqrt(2)*GA_Chroms(minInd, 1)*GA_Chroms(minInd, 2));
            t2 = (M*R)/J;
            T = sqrt(t1^2 + 2*t1*t2*(GA_Chroms(minInd, 2)/(2*R)) + t2^2);
            g1_b = T - Tmax;
            g2_b = Sig - SigMax;
            g3_b = GA_Chroms(minInd, 1) - GA_Chroms(minInd, 4);
            g4_b = 0.125 - GA_Chroms(minInd, 1);
            g5_b = Del - 0.25;
            g6_b = P-Pc;
            g7_b = 0.10471*GA_Chroms(minInd, 1)^2 + 0.04811*GA_Chroms(minInd, 3)*GA_Chroms(minInd, 4)*(14+GA_Chroms(minInd, 2)) - 5;
            
            R = sqrt(((GA_Chroms(maxInd, 2)^2)/4)+((GA_Chroms(maxInd, 1)+GA_Chroms(maxInd, 3))/2)^2);
            Sig = (6*P*L)/(GA_Chroms(maxInd, 4)*GA_Chroms(maxInd, 3)^2);
            Del = (6*P*L^3)/(E*GA_Chroms(maxInd, 4)*GA_Chroms(maxInd, 3)^3);
            Pc = ((4.013*E*sqrt(((GA_Chroms(maxInd, 3)^2)*(GA_Chroms(maxInd, 4)^6))/36))/L^2)*(1-((GA_Chroms(maxInd, 3)/(2*L))*sqrt(E/(4*G))));
            M = P*(L + (GA_Chroms(maxInd, 2)/2));
            J = 2*(sqrt(2)*GA_Chroms(maxInd, 1)*GA_Chroms(maxInd, 2)*(((GA_Chroms(maxInd, 2)^2)/4)+((GA_Chroms(maxInd, 1)+GA_Chroms(maxInd, 3))/2)^2));
            t1 = P/(sqrt(2)*GA_Chroms(maxInd, 1)*GA_Chroms(maxInd, 2));
            t2 = (M*R)/J;
            T = sqrt(t1^2 + 2*t1*t2*(GA_Chroms(maxInd, 2)/(2*R)) + t2^2);
            g1_w = T - Tmax;
            g2_w = Sig - SigMax;
            g3_w = GA_Chroms(maxInd, 1) - GA_Chroms(maxInd, 4);
            g4_w = 0.125 - GA_Chroms(maxInd, 1);
            g5_w = Del - 0.25;
            g6_w = P-Pc;
            g7_w = 0.10471*GA_Chroms(maxInd, 1)^2 + 0.04811*GA_Chroms(maxInd, 3)*GA_Chroms(maxInd, 4)*(14+GA_Chroms(maxInd, 2)) - 5;
            
            mPos = mean(GA_FitVal);
            X = sprintf('Population Converged!\nNumber of Iterations: %d\nBest Value: %0.15f\nBest Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nG5: %0.15f\nG6: %0.15f\nG7: %0.15f\nWorst Value: %0.15f\nWorst Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nG5: %0.15f\nG6: %0.15f\nG7: %0.15f\nMean: %0.15e\n',GA_Curr,GA_FitVal(minInd),GA_Chroms(minInd, :),g1_b,g2_b,g3_b,g4_b,g5_b,g6_b,g7_b,GA_FitVal(maxInd),GA_Chroms(maxInd, :),g1_w,g2_w,g3_w,g4_w,g5_w,g6_w,g7_w,mPos);
            disp(X);
            convRuns = convRuns + 1;
            break;
        end
        PrevDiff = CurrDiff;
    end
    
    if(GA_Curr == GA_MaxItr)
        % if max gen reached
        break;
    end
    
    
    TransPop = zeros(GA_PS, ProbDim);
    
    Arranged = sort(GA_FitVal);
    % Keep half
    for RowNum = 1:ceil(GA_PS/2)
        for i = 1:GA_PS
            if Arranged(RowNum) == GA_FitVal(i)
                TransPop(RowNum, :) = GA_Chroms(i, :);
            end
        end
    end
    
    % Creatwheel
    GA_RouWheel = GA_CreateWheel(GA_PS, GA_FitVal);
    for i = ceil(GA_PS/2)+1:GA_PS
        [Parent1, Parent2] = GA_Selection(GA_PS, GA_Chroms, ProbDim, GA_FitVal, GA_RouWheel);
        SibRep = GA_CrossOver(Parent1, Parent2, GA_cross, ProbDim);
        if rand() <= GA_mut
            SibRep = GA_Mutation(SibRep, DimMinMax, ProbDim);
        end
        TransPop(i, :) = SibRep;
    end
    
    % Make transfer
    GA_Chroms = TransPop;
    GA_Curr = GA_Curr + 1;
end

if GA_Curr >= GA_MaxItr
    for i = 1:GA_PS
        if min(GA_FitVal) == GA_FitVal(i)
            minInd = i;
        end
        if max(GA_FitVal) == GA_FitVal(i)
            maxInd = i;
        end
    end
    Tmax = 13600;
    SigMax = 30000;
    L = 14;
    P = 6000;
    E = 30e6;
    G = 12e6;
    R = sqrt(((GA_Chroms(minInd, 2)^2)/4)+((GA_Chroms(minInd, 1)+GA_Chroms(minInd, 3))/2)^2);
    Sig = (6*P*L)/(GA_Chroms(minInd, 4)*GA_Chroms(minInd, 3)^2);
    Del = (6*P*L^3)/(E*GA_Chroms(minInd, 4)*GA_Chroms(minInd, 3)^3);
    Pc = ((4.013*E*sqrt(((GA_Chroms(minInd, 3)^2)*(GA_Chroms(minInd, 4)^6))/36))/L^2)*(1-((GA_Chroms(minInd, 3)/(2*L))*sqrt(E/(4*G))));
    M = P*(L + (GA_Chroms(minInd, 2)/2));
    J = 2*(sqrt(2)*GA_Chroms(minInd, 1)*GA_Chroms(minInd, 2)*(((GA_Chroms(minInd, 2)^2)/4)+((GA_Chroms(minInd, 1)+GA_Chroms(minInd, 3))/2)^2));
    t1 = P/(sqrt(2)*GA_Chroms(minInd, 1)*GA_Chroms(minInd, 2));
    t2 = (M*R)/J;
    T = sqrt(t1^2 + 2*t1*t2*(GA_Chroms(minInd, 2)/(2*R)) + t2^2);
    g1_b = T - Tmax;
    g2_b = Sig - SigMax;
    g3_b = GA_Chroms(minInd, 1) - GA_Chroms(minInd, 4);
    g4_b = 0.125 - GA_Chroms(minInd, 1);
    g5_b = Del - 0.25;
    g6_b = P-Pc;
    g7_b = 0.10471*GA_Chroms(minInd, 1)^2 + 0.04811*GA_Chroms(minInd, 3)*GA_Chroms(minInd, 4)*(14+GA_Chroms(minInd, 2)) - 5;

    R = sqrt(((GA_Chroms(maxInd, 2)^2)/4)+((GA_Chroms(maxInd, 1)+GA_Chroms(maxInd, 3))/2)^2);
    Sig = (6*P*L)/(GA_Chroms(maxInd, 4)*GA_Chroms(maxInd, 3)^2);
    Del = (6*P*L^3)/(E*GA_Chroms(maxInd, 4)*GA_Chroms(maxInd, 3)^3);
    Pc = ((4.013*E*sqrt(((GA_Chroms(maxInd, 3)^2)*(GA_Chroms(maxInd, 4)^6))/36))/L^2)*(1-((GA_Chroms(maxInd, 3)/(2*L))*sqrt(E/(4*G))));
    M = P*(L + (GA_Chroms(maxInd, 2)/2));
    J = 2*(sqrt(2)*GA_Chroms(maxInd, 1)*GA_Chroms(maxInd, 2)*(((GA_Chroms(maxInd, 2)^2)/4)+((GA_Chroms(maxInd, 1)+GA_Chroms(maxInd, 3))/2)^2));
    t1 = P/(sqrt(2)*GA_Chroms(maxInd, 1)*GA_Chroms(maxInd, 2));
    t2 = (M*R)/J;
    T = sqrt(t1^2 + 2*t1*t2*(GA_Chroms(maxInd, 2)/(2*R)) + t2^2);
    g1_w = T - Tmax;
    g2_w = Sig - SigMax;
    g3_w = GA_Chroms(maxInd, 1) - GA_Chroms(maxInd, 4);
    g4_w = 0.125 - GA_Chroms(maxInd, 1);
    g5_w = Del - 0.25;
    g6_w = P-Pc;
    g7_w = 0.10471*GA_Chroms(maxInd, 1)^2 + 0.04811*GA_Chroms(maxInd, 3)*GA_Chroms(maxInd, 4)*(14+GA_Chroms(maxInd, 2)) - 5;

    mPos = mean(GA_FitVal);
    X = sprintf('Did Not Converge!\nNumber of Iterations: %d\nBest Value: %0.15f\nBest Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nG5: %0.15f\nG6: %0.15f\nG7: %0.15f\nWorst Value: %0.15f\nWorst Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nG5: %0.15f\nG6: %0.15f\nG7: %0.15f\nMean: %0.15e\n',GA_Curr,GA_FitVal(minInd),GA_Chroms(minInd, :),g1_b,g2_b,g3_b,g4_b,g5_b,g6_b,g7_b,GA_FitVal(maxInd),GA_Chroms(maxInd, :),g1_w,g2_w,g3_w,g4_w,g5_w,g6_w,g7_w,mPos);
    disp(X);  
end

    %movie(M,1,120);
    
    if GA_Curr >= GA_MaxItr
        Ans(trials,:) = [GA_Chroms(minInd, :) GA_FitVal(minInd) g1_b g2_b g3_b g4_b g5_b g6_b g7_b 0];
    else % Converged
        Ans(trials,:) = [GA_Chroms(minInd, :) GA_FitVal(minInd) g1_b g2_b g3_b g4_b g5_b g6_b g7_b 1];
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

X = sprintf('\n\nBest OverAll Value: %0.15f\nPosition: %0.15f %0.15f %0.15f %0.15f\nConstraints:\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nG5: %0.15f\nG6: %0.15f\nG7: %0.15f\nMean: %0.15f\nMedian: %0.15f\nStandard Deviation:%0.15f\nWorst Best Overall Value: %0.15f\nNumber of Converged Runs: %d\nRatio of Convergence: %0.15e\nTotal Running Time for all trials: %0.15e\nAverage running time: %0.15e\n', Best, Ans(BesInd, 1), Ans(BesInd, 2), Ans(BesInd, 3), Ans(BesInd, 4), Ans(BesInd, 6), Ans(BesInd, 7), Ans(BesInd, 8), Ans(BesInd, 9),  Ans(BesInd, 10), Ans(BesInd, 11), Ans(BesInd, 12), Mean, Median, StdDev, Worst, convRuns, ConvRatio, totalTime, aveTime);
disp(X);