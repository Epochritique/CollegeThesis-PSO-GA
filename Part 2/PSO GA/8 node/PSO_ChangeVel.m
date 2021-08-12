function [TransPos, TransVel]= PSO_ChangeVel(PopNum, PosPop, VelPop, ProbDim, Pbest, Gbest, w, c1, c2, DimMinMax)
    TransPos = zeros(PopNum, ProbDim);
    TransVel = zeros(PopNum, ProbDim);

    for RowNum = 1:PopNum
        % Calculate new velocity
        TempVel = VelPop(RowNum, 1:ProbDim)*w + c1*rand()*(Pbest(RowNum, 1:ProbDim) - PosPop(RowNum, 1:ProbDim)) + c2*rand()*(Gbest(1, 1:ProbDim) - PosPop(RowNum, 1:ProbDim));
        
        % Check if within Vmin, Vmax
        for i = 1:ProbDim
            VelBound = (DimMinMax(i, 2) - DimMinMax(i, 1))/2;
            if TempVel(i)<-VelBound
                 TempVel(i) = -VelBound;
            elseif TempVel(i)>VelBound
                 TempVel(i) = VelBound;
            end
        end

        TempPos = PosPop(RowNum, :) + TempVel;
        
        uniq=0;
        % Check for Uniqueness
        if length(TempPos(:)) == length(unique(TempPos(:)))
            uniq=1;
        end
        TempPos = Resolve_Uniqueness(TempPos, ProbDim, uniq);
        
                
%         % Check feasibility
%         infcnt=0;
%         for i = 1:ProbDim
%             if TempPos(i)<DimMinMax(i,1) || TempPos(i)>DimMinMax(i,2)
%                 infcnt=infcnt+1;
%             end
%         end
%         
%         % Infeasible
%         if infcnt>0
%             Tao(RowNum, 1) = Tao(RowNum, 1)+1;% Age + 1
%             ConsFac = 2/(abs(-2-Tao(RowNum, 1)-sqrt(Tao(RowNum, 1)^(2) + Tao(RowNum, 1))));
%             TempPos = PosPop(RowNum, :) + ConsFac*TempVel;
%         % Feasible
%         else
%             Tao(RowNum, 1) = 0;
%         end
        
        % Save in respective arrays
        TransPos(RowNum, :) = TempPos;
        TransVel(RowNum, :) = TempVel;
    end
end
