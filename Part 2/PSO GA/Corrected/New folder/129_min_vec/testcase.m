format long;
p = perms([0.1 0.2 0.3 0.4 0.5 0.6 0.7]);

FitVal = [];
for i = 1:120
    FitVal = [FitVal; PSO_GA_Printer(p(i, :), 7, )];
end
ni = min(FitVal);
disp(ni);
x=0;
for i = 1:120
    if abs(FitVal(i)-min(FitVal))<1e-6
        disp(i);
        x = x+1;
    end
end