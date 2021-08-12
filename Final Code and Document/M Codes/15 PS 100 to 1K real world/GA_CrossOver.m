% Function that does the GA crossover mechanic
function SibRep = GA_CrossOver(Parent1, Parent2, GA_cross, ProbDim, DimMinMax, CustCnt)
    CrossVal = rand();
    % Cross or Not
    if CrossVal <= GA_cross
        % Cross at random point
        CrossPos = randi(ProbDim-1);
        Sibling1 = [Parent1(1:CrossPos) Parent2(CrossPos+1:ProbDim)];
        Sibling2 = [Parent2(1:CrossPos) Parent1(CrossPos+1:ProbDim)];
    else
        % No crossing
        Sibling1 = Parent1;
        Sibling2 = Parent2;
    end
%     uniq=0;
%     % Check for Uniqueness
%     if length(Sibling1(:)) == length(unique(Sibling1(:)))
%         uniq=1;
%     end
%     Sibling1 = Resolve_Uniqueness(Sibling1, ProbDim, uniq);

%     uniq=0;
%     % Check for Uniqueness
%     if length(Sibling2(:)) == length(unique(Sibling2(:)))
%         uniq=1;
%     end
%     Sibling2 = Resolve_Uniqueness(Sibling2, ProbDim, uniq);

    
    % Compute Fitness Values
    [Sib1_FitVal,Sib1_Exc]= PSO_GA_Eval(Sibling1, ProbDim, CustCnt);
    [Sib2_FitVal,Sib2_Exc]= PSO_GA_Eval(Sibling2, ProbDim, CustCnt);    
    if Sib1_FitVal <= Sib2_FitVal
        if Sib2_Exc < Sib1_Exc
            SibRep = Sibling2;
        else
            SibRep = Sibling1;
        end
    else
        if Sib1_Exc < Sib2_Exc
            SibRep = Sibling1;
        else
            SibRep = Sibling2;
        end
    end
end
