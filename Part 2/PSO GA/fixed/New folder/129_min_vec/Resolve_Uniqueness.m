% Function for checking and resolving uniqueness of values for each
% candidate solution
function TempPos = Resolve_Uniqueness(TempPos, ProbDim, uniq)
    while uniq==0
        % Resolve
        for i = 1:ProbDim
            % Check for Uniqueness
            if ~(length(TempPos(1:i)) == length(unique(TempPos(1:i))))
                % Adjust values using small changes
                if rand()>0.5
                    TempPos(i) = TempPos(i)+rand()*0.1;
                else
                    TempPos(i) = TempPos(i)-rand()*0.1;
                end
            end
        end
        % Check for Uniqueness
        if length(TempPos(:)) == length(unique(TempPos(:)))
            uniq=1;
        end
    end
end