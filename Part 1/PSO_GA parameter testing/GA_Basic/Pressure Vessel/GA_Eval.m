function FitVal = GA_Eval(Chrom, ProbDim)
    FitVal = 0.6224*Chrom(1)*Chrom(3)*Chrom(4) + 1.7781*Chrom(2)*(Chrom(3)^2) + 3.1661*(Chrom(1)^2)*Chrom(4) + 19.84*(Chrom(1)^2)*Chrom(3);
end
