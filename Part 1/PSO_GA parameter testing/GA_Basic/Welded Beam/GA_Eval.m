function FitVal = GA_Eval(Chrom, ProbDim)
    FitVal = 1.10471*(Chrom(1)^2)*Chrom(2) + 0.04811*Chrom(3)*Chrom(4)*(14+Chrom(2));
end
