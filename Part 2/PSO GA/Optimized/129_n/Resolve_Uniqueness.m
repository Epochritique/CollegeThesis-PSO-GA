function TempPos = Resolve_Uniqueness(TempPos, ProbDim, uniq)
    while uniq==0
        % Resolve
        for i = 1:ProbDim
            % Check for Uniqueness
            if ~(length(TempPos(1:i)) == length(unique(TempPos(1:i))))
                if rand()>0.5
                    TempPos(i) = TempPos(i)+0.0001;
                else
                    TempPos(i) = TempPos(i)-0.0001;
                end
            end
        end
        % Check for Uniqueness
        if length(TempPos(:)) == length(unique(TempPos(:)))
            uniq=1;
        end
    end
end