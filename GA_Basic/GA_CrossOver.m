function SibRep = GA_CrossOver(Parent1, Parent2, GA_cross, ProbDim)
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
    
    % Compute Fitness Values
    Sib1_FitVal = GA_Eval(Sibling1, ProbDim);
    Sib2_FitVal = GA_Eval(Sibling2, ProbDim);
    if Sib1_FitVal >= Sib2_FitVal
        SibRep = Sibling1;
    else
        SibRep = Sibling2;
    end
end
