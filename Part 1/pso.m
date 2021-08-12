%LICNACHAN, LANCE OLIVER C.
%2014-64880

format long;
PopNum = 50;
Bound_upper = 100;
Bound_lower = -100;
Bound_upper_velo = 3;
Bound_lower_velo = -3;
MaxItr = 100000;
Threshold = 500;
c1 = 2;
c2 = 2;
w = 0.9;

InitPop = zeros(PopNum, 2);
TransPop = zeros(PopNum, 2);
RowNum = 1;
FitVal = zeros(PopNum, 1);
Pbest = zeros(PopNum, 2);
Gbest = zeros(1,2);


% Get initial population
for RowNum = 1:PopNum
    TempPos = Bound_upper + (2*(Bound_lower))*rand(1,1,'double'); %Randomize a position
    TempVel = Bound_upper_velo + (2*(Bound_lower_velo))*rand(1,1,'double'); %Randomize a velocity
    InitPop(RowNum, :) = [TempPos, TempVel];  %Add to population
end

% Evaluate each particle's value
for RowNum = 1:PopNum
    FitVal(RowNum, 1) = ObjFxn(InitPop(RowNum,1));
end

% Get Initial Pbest and Gbest
for RowNum = 1:PopNum
    Pbest(RowNum, 1) = InitPop(RowNum,1);
    Pbest(RowNum, 2) = FitVal(RowNum,1);
end

Gbest(2) = min(FitVal);
for RowNum = 1:PopNum
    if Pbest(RowNum,2) == Gbest(2)
        Gbest(1) = Pbest(RowNum, 1);
        break;
    end
end

CurrItr=1;
while CurrItr <= MaxItr
    if(CurrItr == Threshold)
        w = 0.4;
    end
    % Evaluate each particle's value
    for RowNum = 1:PopNum
        FitVal(RowNum, 1) = ObjFxn(InitPop(RowNum,1));
    end

    % Get Pbest and Gbest
    for RowNum = 1:PopNum
        if Pbest(RowNum, 2) > FitVal(RowNum,1)
            Pbest(RowNum, 1) = InitPop(RowNum,1);
            Pbest(RowNum, 2) = FitVal(RowNum,1);
        end
    end
    if Gbest(2) > min(FitVal)
        Gbest(2) = min(FitVal);
        for RowNum = 1:PopNum
            if Pbest(RowNum,2) == Gbest(2)
                Gbest(1) = Pbest(RowNum, 1);
                break;
            end
        end
    end
    
    % Change Velocity According to Best Positions
    for RowNum = 1:PopNum
        TempVel = InitPop(RowNum, 2)*w + c1*rand()*(Pbest(RowNum, 1) - InitPop(RowNum, 1)) + c2*rand()*(Gbest(1) - InitPop(RowNum, 1));
        TempPos = InitPop(RowNum, 1) + TempVel;
        TransPop(RowNum, :) = [TempPos, TempVel];  %Add to population
    end

    ctr = 0;
    for RowNum = 1:PopNum
        if Gbest(2) == FitVal(RowNum,1)
            ctr = ctr+1;
        end
    end
        
    if ctr == PopNum && (Gbest(2) < 1e-7)
        X = sprintf('Converged!\nNumber of Iterations: %d\nBest Value: %d\nBest Location: %d',CurrItr,Gbest(2),Gbest(1));
        disp(X);
        break;
    end
    InitPop = TransPop;
    CurrItr=CurrItr+1;
end

if CurrItr > MaxItr
    X = sprintf('Did not Converge!\nNumber of Iterations: %d\nBest Value: %d\nBest Location: %d\nNumber of Converged Individuals: %d',CurrItr,Gbest(2),Gbest(1),ctr);
    disp(X);
end