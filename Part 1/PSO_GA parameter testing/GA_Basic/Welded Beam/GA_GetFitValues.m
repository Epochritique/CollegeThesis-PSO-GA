function GA_FitVal = GA_GetFitValues(GA_PS, GA_Chroms, ProbDim, DimMinMax)
    GA_FitVal = zeros(GA_PS, 1);
    
    % Evaluate each particle's value
    for RowNum = 1:GA_PS
        GA_FitVal(RowNum, 1) = GA_Eval(GA_Chroms(RowNum, :), ProbDim);
    end
    
    FitWorst = max(GA_FitVal);
    for RowNum = 1:GA_PS
        probs = 0;
        
        for i = 1:ProbDim
           if ~ (GA_Chroms(RowNum, i) >= DimMinMax(i, 1) && GA_Chroms(RowNum, i) <= DimMinMax(i, 2))
               probs=probs+1;% compute the number of boundary errors
           end
        end

    Tmax = 13600;
    SigMax = 30000;
    L = 14;
    P = 6000;
    E = 30e6;
    G = 12e6;
    R = sqrt(((GA_Chroms(RowNum, 2)^2)/4)+((GA_Chroms(RowNum, 1)+GA_Chroms(RowNum, 3))/2)^2);
    Sig = (6*P*L)/(GA_Chroms(RowNum, 4)*GA_Chroms(RowNum, 3)^2);
    Del = (6*P*L^3)/(E*GA_Chroms(RowNum, 4)*GA_Chroms(RowNum, 3)^3);
    Pc = ((4.013*E*sqrt(((GA_Chroms(RowNum, 3)^2)*(GA_Chroms(RowNum, 4)^6))/36))/L^2)*(1-((GA_Chroms(RowNum, 3)/(2*L))*sqrt(E/(4*G))));
    M = P*(L + (GA_Chroms(RowNum, 2)/2));
    J = 2*(sqrt(2)*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 2)*(((GA_Chroms(RowNum, 2)^2)/4)+((GA_Chroms(RowNum, 1)+GA_Chroms(RowNum, 3))/2)^2));
    t1 = P/(sqrt(2)*GA_Chroms(RowNum, 1)*GA_Chroms(RowNum, 2));
    t2 = (M*R)/J;
    T = sqrt(t1^2 + 2*t1*t2*(GA_Chroms(RowNum, 2)/(2*R)) + t2^2);
    % Check if it satisfies constraint equations
    g1 = T - Tmax;
    g2 = Sig - SigMax;
    g3 = GA_Chroms(RowNum, 1) - GA_Chroms(RowNum, 4);
    g4 = 0.125 - GA_Chroms(RowNum, 1);
    g5 = Del - 0.25;
    g6 = P-Pc;
    g7 = 0.10471*GA_Chroms(RowNum, 1)^2 + 0.04811*GA_Chroms(RowNum, 3)*GA_Chroms(RowNum, 4)*(14+GA_Chroms(RowNum, 2)) - 5;
        if ~(g1<=0)
            probs=probs+1;% compute the number of errors
        end
        if ~(g2<=0)    
            probs=probs+1;% compute the number of errors
        end
        if ~(g3<=0)   
            probs=probs+1;% compute the number of errors
        end
        if ~(g4<=0)   
            probs=probs+1;% compute the number of errors
        end
        if ~(g5<=0)    
            probs=probs+1;% compute the number of errors
        end
        if ~(g6<=0)    
            probs=probs+1;% compute the number of errors
        end
        if ~(g7<=0)    
            probs=probs+1;% compute the number of errors
        end
        if(probs~=0)
            GA_FitVal(RowNum, 1) = FitWorst + abs(g1) + abs(g2) + abs(g3) + abs(g4) + abs(g5) + abs(g6) + abs(g7);
        end
    end
end
