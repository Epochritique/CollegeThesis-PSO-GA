%   LICNACHAN, LANCE OLIVER C.
%   2014-64880
%   June 2018

format long;
CurPopNum = 100;
maxitrun=1;
while maxitrun <= 5
    Y = sprintf('\nCurPopNum: %d\n',CurPopNum);
    disp(Y);
    Y = sprintf('\nMaxIt: %d\n',(CurPopNum*(10*maxitrun)));
    disp(Y);
    
    CustCnt = 129;                  % Number of Collection sites
    ProbDim = (2*CustCnt)-1;        % Problem Dimensions
    RunMax = 10;                    % Number of Trials
    convRuns=0;
    Ans = zeros(RunMax, ProbDim++1+1);
    timeRec = zeros(RunMax,1);

    for trials = 1:RunMax
    tic;    % Start Timer
    Y = sprintf('\nTrial: %d\n',trials);
    disp(Y);

    % Variables specific to the problem
    DimMinMax = zeros(ProbDim, 2);
    for i=1:ProbDim
        DimMinMax(i, :) = [0 1];
    end

    % Variables specific to the algorithm
    AcceptThreshold = 1e-5;
    
    % PSO Variables
    PopNum = CurPopNum;
    PSO_Curr = 1;
    PSO_Max = CurPopNum*10*maxitrun;
    c1 = 1.5;
    c2 = 1.5;
    wmax = 0.9;
    wmin = 0.4;
    
    % GA Variables
    GA_cross = 0.85;
    GA_mut = 0.02;
    GA_y = 10;
    GA_B = 15;
    GA_NumMax = ceil(PopNum*0.1);
    GA_NumMin = 1;
    GA_MinPS = 20;
    GA_MaxPS = 10;
    GA_MinItr = 10;
    GA_MaxItr = 15;
    GA_MaxItr = floor(GA_MinItr + ((PSO_Curr/PSO_Max)^GA_B)*(GA_MaxItr-GA_MinItr));
    GA_Num = floor(GA_NumMax - ((PSO_Curr/PSO_Max)^GA_y)*(GA_NumMax-GA_NumMin));
    GA_PS = floor(GA_MinPS + ((PSO_Curr/PSO_Max)^GA_y)*(GA_MaxPS-GA_MinPS));

    TransPos = zeros(PopNum, ProbDim);
    TransVel = zeros(PopNum, ProbDim);
    % Initialization Step
    [PosPop, VelPop] = PSO_InitPop(PopNum, ProbDim, DimMinMax);
    FitVal = PSO_GA_GetFitValues(PopNum, PosPop, ProbDim, CustCnt);
    [Pbest, Gbest]= PSO_GetInitPGBest(PopNum, PosPop, ProbDim, FitVal);

    PrevDiff = 0;
    while PSO_Curr <= PSO_Max
        % Evaluate Fitness
        FitVal = PSO_GA_GetFitValues(PopNum, PosPop, ProbDim, CustCnt);

        % Check if termination conditions have been met
        if(PSO_Curr == 1)
            PrevDiff = max(FitVal) - min(FitVal);
        else
            CurrDiff = max(FitVal) - min(FitVal);
            % Check for population convergence
            if PrevDiff - CurrDiff <= AcceptThreshold && CurrDiff <= AcceptThreshold && PSO_Curr < PSO_Max
                for i = 1:PopNum
                    if min(FitVal) == FitVal(i)
                        minInd = i;
                    end
                    if max(FitVal) == FitVal(i)
                        maxInd = i;
                    end
                end
                mPos = mean(FitVal);
                seq_route = ObtainSequence(PosPop(minInd,:), ProbDim, CustCnt);
                disp(seq_route);
                X = sprintf('Population Converged!\nNumber of Iterations: %d\nBest Value: %0.15f\nWorst Value: %0.15f\nMean: %0.15f',PSO_Curr,FitVal(minInd),FitVal(maxInd),mPos);
                disp(X);
                convRuns = convRuns + 1;
                break;
            end
            PrevDiff = CurrDiff;
        end

        if(PSO_Curr == PSO_Max)
            % if max gen reached
            break;
        end

        % Get best values
        [Pbest, Gbest] = PSO_GetPGBest(PopNum, PosPop,  ProbDim, FitVal, Pbest, Gbest);

        % Change value according to how current iteration
        w = wmax-(wmax-wmin)*(PSO_Curr/PSO_Max);

        % Calculate new velocities and move
        [TransPos, TransVel] = PSO_ChangeVel(PopNum, PosPop, VelPop, ProbDim, Pbest, Gbest, w, c1, c2, DimMinMax);

        % Evaluate
        TransFitVal = PSO_GA_GetFitValues(PopNum, TransPos, ProbDim, CustCnt);

        % GA Portion
        PSO_Arranged = sort(TransFitVal);
        GA_Num_Curr = 1;
        while GA_Num_Curr <= GA_Num
            % Get one from best individuals
            for RowNum = 1:PopNum
                if TransFitVal(RowNum) == PSO_Arranged(GA_Num_Curr)
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
                GA_FitVal = PSO_GA_GetFitValues(GA_PS, GA_Chroms, ProbDim, CustCnt);
                TransPop = zeros(GA_PS, ProbDim);

                % Keep Elite
                % Note: we keep the elite feasible
                for i = 1:GA_PS
                    [~, TotExc1] = PSO_GA_Eval(GA_Fit_Chrom, ProbDim, (ProbDim+1)/2);
                    [~, TotExc2] = PSO_GA_Eval(GA_Chroms(i,:), ProbDim, (ProbDim+1)/2);
                    if GA_Fit_Elite > GA_FitVal(i) && TotExc1 >= TotExc2
                        GA_Fit_Chrom = GA_Chroms(i,:);
                        GA_Fit_Elite = GA_FitVal(i);
                    end
                end

                % Exit at max iterations
                if GA_Curr == GA_MaxItr
                    break;
                end

                % Create Wheel
                GA_RouWheel = GA_CreateWheel(GA_PS, GA_FitVal);

                % Create the next population
                for i = 1:GA_PS
                    % Select 2 Parents
                    [Parent1, Parent2] = GA_Selection(GA_PS, GA_Chroms, ProbDim, GA_FitVal, GA_RouWheel);
                    % Cross-over
                    SibRep = GA_CrossOver(Parent1, Parent2, GA_cross, ProbDim, DimMinMax, CustCnt);
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
            % Replace the infeasible/worst individuals
            for RowNum = 1:PopNum
                if TransFitVal(RowNum) == PSO_Arranged(PopNum-GA_Num_Curr+1);
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

        % Update PSO variables
        PosPop = TransPos;
        VelPop = TransVel;
        PSO_Curr = PSO_Curr + 1;
    end

    % Some printing
    if PSO_Curr >= PSO_Max
        for i = 1:PopNum
            if min(FitVal) == FitVal(i)
                minInd = i;
            end
            if max(FitVal) == FitVal(i)
                maxInd = i;
            end
        end
        seq_route = ObtainSequence(PosPop(minInd,:), ProbDim, CustCnt);
        disp(seq_route);
        mPos = mean(FitVal);
        X = sprintf('Did Not Converge!\nNumber of Iterations: %d\nBest Value: %0.15f\nWorst Value: %0.15f\nMean: %0.15e\n',PSO_Curr,FitVal(minInd),FitVal(maxInd),mPos);
        disp(X);  
    end

        if PSO_Curr >= PSO_Max
            Ans(trials,:) = [PosPop(minInd, :) FitVal(minInd) 0];
        else % Converged
            Ans(trials,:) = [PosPop(minInd, :) FitVal(minInd) 1];
        end

        timeRec(trials) = toc;
        X = sprintf('Running Time for this trial: %0.15e\n', timeRec(trials));
        disp(X);
    end
    if convRuns > 0
        % Get Best Fit
        ConvVals = zeros(convRuns,1);
        i=1;
        for o = 1:RunMax
            if(Ans(o,ProbDim++1+1) == 1)
                ConvVals(i) = Ans(o,ProbDim+1);
                i=i+1;
            end
        end

        Best = min(ConvVals);
        for o = 1:RunMax
            if min(ConvVals) == Ans(o,ProbDim+1)
                BesInd = o;
                break;
            end
        end

        % Generate Stats
        Mean = mean(ConvVals);
        StdDev = std(ConvVals);
        Median = median(ConvVals);
        Worst = max(ConvVals);

    else
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
        for o = 1:RunMax
            if min(Vals) == Ans(o,ProbDim+1)
                BesInd = o;
            end
        end
    end

    ConvRatio = convRuns/RunMax;
    totalTime = sum(timeRec);
    aveTime = mean(timeRec);

    X = sprintf('\n\nBest OverAll Value: %0.15f\nMean: %0.15f\nMedian: %0.15f\nStandard Deviation:%0.15f\nWorst Best Overall Value: %0.15f\nNumber of Converged Runs: %d\nRatio of Convergence: %0.15e\nTotal Running Time for all trials: %0.15e\nAverage running time: %0.15e\n', Best, Mean, Median, StdDev, Worst, convRuns, ConvRatio, totalTime, aveTime);
    disp(X);

    save(sprintf('run_%d_%d_%d',CustCnt,PopNum, PSO_Max), 'Ans', 'Best', 'Mean', 'Median', 'StdDev', 'Worst', 'convRuns', 'ConvRatio', 'totalTime', 'aveTime');

    maxitrun=maxitrun+1;
end