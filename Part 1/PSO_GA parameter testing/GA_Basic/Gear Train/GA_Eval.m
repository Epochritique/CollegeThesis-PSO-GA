function FitVal = GA_Eval(Chrom, ProbDim)
    FitVal = ((1/6.931) - ((Chrom(1)*Chrom(2))/(Chrom(3)*Chrom(4))))^2;
end
