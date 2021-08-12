%LICNACHAN, LANCE OLIVER C.
%2014-64880

format long;
% rng('shuffle');
rune=1;
while rune <= 3
    CurPopNum = 150;
    if rune == 2
        CurPopNum = 250;
    elseif rune == 3
        CurPopNum = 500;
    end
    Y = sprintf('\nCurPopNum: %d\n',CurPopNum);
    disp(Y);
    
    mun=1;
    while mun <= 5
    Y = sprintf('\nMaxIt: %d\n',(CurPopNum*(2*mun)));
    disp(Y);
    
        ProbDim = 4;
        ConsNum = 4;
        convRuns = 0;
        RunMax = 10;
        Ans = zeros(RunMax, ProbDim+ConsNum+1+1);
        timeRec = zeros(RunMax,1);
        for trials = 1:RunMax
        tic;
        Y = sprintf('\nTrial: %d\n',trials);
        disp(Y);

        % Variables specific to the problem

        DimMinMax = zeros(ProbDim, 2);
        DimMinMax(1, :) = [1 99]*0.0625;
        DimMinMax(2, :) = [1 99]*0.0625;
        DimMinMax(3, :) = [10 200];
        DimMinMax(4, :) = [10 200];

        % Variables specific to the algorithm
        AcceptThreshold = 1e-5;
        PopNum = CurPopNum;
        % PopNum = 200;
        PSO_Curr = 1;
        PSO_Max = CurPopNum*(2*mun);
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

            if(PSO_Curr == 1)
                PrevDiff = max(FitVal) - min(FitVal);
            else
                CurrDiff = max(FitVal) - min(FitVal);
                % disp(CurrDiff);
                % Check for population convergence
                if PrevDiff - CurrDiff < AcceptThreshold && CurrDiff < AcceptThreshold && PSO_Curr < PSO_Max
                    for i = 1:PopNum
                        if min(FitVal) == FitVal(i)
                            minInd = i;
                        end
                        if max(FitVal) == FitVal(i)
                            maxInd = i;
                        end
                    end
                    g1_b = -PosPop(minInd, 1) + 0.0193*PosPop(minInd, 3);
                    g2_b = -PosPop(minInd, 2) + 0.0095*PosPop(minInd, 3);
                    g3_b = -(pi*(PosPop(minInd, 3)^2)*PosPop(minInd, 4)) - ((4/3)*pi*(PosPop(minInd, 3)^3)) + 1296000;
                    g4_b = PosPop(minInd, 4) - 240;

                    g1_w = -PosPop(maxInd, 1) + 0.0193*PosPop(maxInd, 3);
                    g2_w = -PosPop(maxInd, 2) + 0.0095*PosPop(maxInd, 3);
                    g3_w = -(pi*(PosPop(maxInd, 3)^2)*PosPop(maxInd, 4)) - ((4/3)*pi*(PosPop(maxInd, 3)^3)) + 1296000;
                    g4_w = PosPop(maxInd, 4) - 240;

                    mPos = mean(FitVal);
                    X = sprintf('Population Converged!\nNumber of Iterations: %d\nBest Value: %0.15f\nBest Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nWorst Value: %0.15f\nWorst Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nMean: %0.15e\n',PSO_Curr,FitVal(minInd), PosPop(minInd, :), g1_b,g2_b,g3_b,g4_b, FitVal(maxInd), PosPop(maxInd, :), g1_w,g2_w,g3_w,g4_w, mPos);
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
            [Pbest, Gbest] = PSO_GetPGBest(PopNum, PosPop, ProbDim, FitVal, Pbest, Gbest);

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
                    for i = 1:GA_PS
                        if GA_Fit_Elite < GA_FitVal(i)
                            GA_Chroms(i,:)=GA_Fit_Chrom;
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
                % Obtain current best
                Arranged = sort(GA_FitVal);
                if Arranged(1) < GA_Fit_Elite
                    GA_Fit_Elite = Arranged(1);
                    for i = 1:GA_PS
                        if Arranged(1) == GA_FitVal(i)
                            GA_Fit_Chrom = GA_Chroms(i,:);
                        end
                    end
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
            g1_b = -PosPop(minInd, 1) + 0.0193*PosPop(minInd, 3);
            g2_b = -PosPop(minInd, 2) + 0.0095*PosPop(minInd, 3);
            g3_b = -(pi*(PosPop(minInd, 3)^2)*PosPop(minInd, 4)) - ((4/3)*pi*(PosPop(minInd, 3)^3)) + 1296000;
            g4_b = PosPop(minInd, 4) - 240;

            g1_w = -PosPop(maxInd, 1) + 0.0193*PosPop(maxInd, 3);
            g2_w = -PosPop(maxInd, 2) + 0.0095*PosPop(maxInd, 3);
            g3_w = -(pi*(PosPop(maxInd, 3)^2)*PosPop(maxInd, 4)) - ((4/3)*pi*(PosPop(maxInd, 3)^3)) + 1296000;
            g4_w = PosPop(maxInd, 4) - 240;

            mPos = mean(FitVal);
            X = sprintf('Did Not Converge!\nNumber of Iterations: %d\nBest Value: %0.15f\nBest Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nWorst Value: %0.15f\nWorst Position: %0.15f %0.15f %0.15f %0.15f\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nMean: %0.15e\n',PSO_Curr,FitVal(minInd), PosPop(minInd, :), g1_b,g2_b,g3_b,g4_b, FitVal(maxInd), PosPop(maxInd, :), g1_w,g2_w,g3_w,g4_w, mPos);
            disp(X);
        end
            if PSO_Curr >= PSO_Max
                Ans(trials,:) = [PosPop(minInd, :) FitVal(minInd) g1_b g2_b g3_b g4_b 0];
            else % Converged
                Ans(trials,:) = [PosPop(minInd, :) FitVal(minInd) g1_b g2_b g3_b g4_b 1];
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
                if(Ans(o,ProbDim+ConsNum+1+1) == 1)
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

        X = sprintf('\n\nBest OverAll Value: %0.15f\nPosition: %0.15f %0.15f %0.15f %0.15f\nConstraints:\nG1: %0.15f\nG2: %0.15f\nG3: %0.15f\nG4: %0.15f\nMean: %0.15f\nMedian: %0.15f\nStandard Deviation:%0.15f\nWorst Best Overall Value: %0.15f\nNumber of Converged Runs: %d\nRatio of Convergence: %0.15e\nTotal Running Time for all trials: %0.15e\nAverage running time: %0.15e\n', Best, Ans(BesInd, 1), Ans(BesInd, 2), Ans(BesInd, 3), Ans(BesInd, 4), Ans(BesInd, 6), Ans(BesInd, 7), Ans(BesInd, 8), Ans(BesInd, 9), Mean, Median, StdDev, Worst, convRuns, ConvRatio, totalTime, aveTime);
        disp(X);

        mun=mun+1;
    end
    rune = rune+1;
end
