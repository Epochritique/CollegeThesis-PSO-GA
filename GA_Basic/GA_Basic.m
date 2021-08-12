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
DimMinMax(2, :) = [0 1];
DimMinMax(3, :) = [0 1];
DimMinMax(4, :) = [0 1];
DimMinMax(5, :) = [0 1];

% Variables specific to the algorithm
AcceptThreshold = 1e-6;
GA_cross = 0.85;
GA_mut = 0.02;
GA_Curr = 1;
GA_MaxItr = 1000;
GA_PS = 20; % Population Size

GA_Chroms = GA_InitPop(GA_PS, ProbDim, DimMinMax);

PrevDiff = 0;
while GA_Curr <= GA_MaxItr
    GA_FitVal = GA_GetFitValues(GA_PS, GA_Chroms, ProbDim);
    TransPop = zeros(GA_PS, ProbDim);
    
    Arranged = sort(GA_FitVal);
    for RowNum = 1:ceil(GA_PS/2)
        for i = 1:GA_PS
            if Arranged(RowNum) == GA_FitVal(i)
                TransPop(RowNum, :) = GA_Chroms(i, :);
            end
        end
    end
    
    GA_RouWheel = GA_CreateWheel(GA_PS, GA_FitVal);
    for i = ceil(GA_PS/2)+1:GA_PS
        [Parent1, Parent2] = GA_Selection(GA_PS, GA_Chroms, ProbDim, GA_FitVal, GA_RouWheel);
        SibRep = GA_CrossOver(Parent1, Parent2, GA_cross, ProbDim);
        if rand() <= GA_mut
            SibRep = GA_Mutation(SibRep, DimMinMax, ProbDim);
        end
        TransPop(i, :) = SibRep;
    end
    GA_Chroms = TransPop;
        
    for i = 1:GA_PS
        if min(GA_FitVal) == GA_FitVal(i)
            BestInd = GA_Chroms(i, :);
        end
    end
    if(GA_Curr == 1)
        PrevDiff = max(GA_FitVal) - min(GA_FitVal);
    else
        CurrDiff = max(GA_FitVal) - min(GA_FitVal);
        if (PrevDiff - CurrDiff < AcceptThreshold) && CurrDiff < AcceptThreshold
            X = sprintf('Population Converged!\nNumber of Iterations: %d\nBest Value: %d\nBest Invidividual: %d %d %d %d %d\n',GA_Curr,min(GA_FitVal),BestInd);
            disp(X);
            break;
        end
        PrevDiff = CurrDiff;
    end

    GA_Curr = GA_Curr + 1;
end

if GA_Curr > GA_MaxItr
    X = sprintf('Did Not Converge!\nNumber of Iterations: %d\nBest Value: %d\n',GA_Curr,min(GA_FitVal));
    disp(X);
end
end