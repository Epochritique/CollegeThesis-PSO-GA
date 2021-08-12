function Chrom = GA_Mutation(Chrom, DimMinMax, ProbDim)
    MutPos = randi(ProbDim);    
    % Non-Binary
    Chrom(1, MutPos) = randi([DimMinMax(MutPos, 1), DimMinMax(MutPos, 2)], 1);
end
