function Chrom = GA_Mutation(Chrom, DimMinMax, ProbDim)
    MutPos = randi(ProbDim);    
    % Non-Binary
    Chrom(1, MutPos) = (DimMinMax(MutPos, 2)- DimMinMax(MutPos, 1))*rand(1,'double') + DimMinMax(MutPos, 1);
    
    uniq=0;
    while uniq==0
        % Check for Uniqueness
        if length(Chrom(:)) == length(unique(Chrom(:)))
            uniq=1;
        else
            Chrom(1, MutPos) = (DimMinMax(MutPos, 2)- DimMinMax(MutPos, 1))*rand(1,'double') + DimMinMax(MutPos, 1);
        end
    end

end
