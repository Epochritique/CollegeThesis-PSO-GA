%LICNACHAN, LANCE OLIVER C.
%2014-64880

%Global constants
tic;
AllNum = 5;             %number of alleles 
PopNum = 50;           %number of population
GenNum = PopNum*10000; 
Bound_up = 10;        %upper bound
Bound_do = -10;       %lower bound
MutRate = .2;           %Mutation Rate
CrosRate = .8;         %Crossover Rate
Reps = 5;

%INITIAL POPULATION
InitPop = zeros(PopNum, AllNum);
RowNum = 1;
for RowNum = 1:PopNum
    Temp = (Bound_up+0.00001) + (2*(Bound_do+0.00001))*rand(1,'double'); %Randomize a chromosome\
    TempRow = zeros(1,AllNum);
    for p = 1:AllNum
        TempRow(1,p) = Temp;
    end
    InitPop(RowNum, 1:AllNum) = TempRow;  %Add to population
end

CurGen = 1;
while CurGen <= GenNum
    %SOLVE FOR OBJECTIVE FITNESS
    ObjFitVal = zeros(PopNum, 1);
    for RowNum = 1:PopNum
        num=0;
        for p = 1:AllNum
            num = num + InitPop(RowNum,p);
        end
        ObjFitVal(RowNum, 1) = ObjFxn(num);
    end
    
    Arranged = sort(ObjFitVal); %Arrange the objective fitness values
    if CurGen == 1
        min = Arranged(1);
        counter = 0;
    end
    MaxRow = 0;
    RowNum = 1;
    while RowNum <= PopNum %Find the top chromosome
        if ObjFitVal(RowNum,1) == Arranged(1)
            MaxRow = RowNum;
            Top1(1:AllNum) = InitPop(RowNum, 1:AllNum);
        end
        RowNum = RowNum + 1;
    end
    if min > Arranged(1)
        min = Arranged(1);
        counter = 0;
    elseif min == Arranged(1)
        counter = counter+1;
    end
    if(5e-5>Arranged(1) && counter == Reps)%When minimum value has been reached
        TimeEl = toc;
        disp('Min Value Found!');
        break;
    end
    
    if CurGen == GenNum %When at the end of N generations
        TimeEl = toc;
        disp('Maximum Generations Reached!');
        break;
    end
   
    %SELECTION
    InterPop = zeros(PopNum, AllNum);
    RowNum = 1;
    %Keep the best 50%
    while RowNum < round(PopNum/2)
        Ctr = 1;
        while Ctr <= PopNum
            if Arranged(RowNum)== ObjFitVal(Ctr)
                InterPop(RowNum, 1:AllNum) = InitPop(Ctr, 1:AllNum);
            end
             Ctr=Ctr+1;
        end
        RowNum = RowNum + 1;
    end
    RowNum = round(PopNum/2);
    %Generate rest by tournament
    while RowNum <= PopNum
        Fighter1 = randi([1 PopNum]); % Generate 1st fighter
        Fighter2 = Fighter1;
        while (Fighter2 == Fighter1) %Generate 2nd fighter that is not the same individual as the 1st fighter
            Fighter2 = randi([1 PopNum]);
        end
        if(ObjFitVal(Fighter1)<ObjFitVal(Fighter2)) %if the random number is less than the probability of the 1st fighter, place it in the intermediate population 
            InterPop(RowNum, 1:AllNum) = InitPop(Fighter1, 1:AllNum);
        else %else place the 2nd fighter in the intermediate population 
            InterPop(RowNum, 1:AllNum) = InitPop(Fighter2, 1:AllNum);
        end
        RowNum = RowNum + 1;
    end
    %CROSSING
    RowNum = 1;
    ParNum = 1;
    while RowNum <= PopNum
        GenRand = rand;     %Generate a random number from 0 to 1 for each chromosome
        if GenRand < CrosRate   %Select the chromosome for which the random number is within the crossing rate
            Parents(ParNum, 1:AllNum) = InterPop(RowNum, 1:AllNum);
            Rows(ParNum) = RowNum;
            ParNum = ParNum + 1;    %Incerement the size of the array of parents
        end
        RowNum = RowNum + 1;%Traverse the new population
    end
    ParNum = ParNum - 1;    %reduce the total increment to match the size of the array
    if ParNum ~= 1  %IF JUST 1 PARENT, skip crossing process
        RowNum = 1;
        while RowNum <= ParNum %Generate a random number to find a mate for each parent with another in the parent pool
            GenRand = randi([1 ParNum]);
            if GenRand ~= RowNum
                ParMat(RowNum, 1:AllNum) = Parents(GenRand, 1:AllNum);
                CrossPnt(RowNum, 1) = randi([1 AllNum-1]); %Generate a random crossing point for each pair
                RowNum = RowNum + 1;
            end
        end
        
        %Proceed with crossing
        RowNum = 1;
        while RowNum <= ParNum
            InterPop(Rows(RowNum), 1:CrossPnt(RowNum)) = Parents(RowNum, 1:CrossPnt(RowNum)); %Set the new alleles of the chromosome
            InterPop(Rows(RowNum), CrossPnt(RowNum)+1:AllNum) = ParMat(RowNum, CrossPnt(RowNum)+1:AllNum);  %Set the new alleles of the chromosome
            RowNum = RowNum + 1;
        end
    end
    
    %Mutation
    while RowNum <= PopNum %Run all the individuals
        Deter = rand; %Determine if individual will undergo Mutation
        if Deter <= MutRate    %If generated number is within mutation rate, mutate
            Bit = randi([1 AllNum]);    %Randomize which bit to mutate
            InterPop(RowNum, Bit) = ((Bound_up+0.00001) + (2*(Bound_do+0.00001))*rand(1,'double'))/5; %Set the new value
        end
        RowNum = RowNum + 1;
    end
    InitPop = InterPop;
    CurGen = CurGen + 1;
end
disp('Generations: ');
disp(CurGen);
disp('Min Val: ');
disp(ObjFitVal(MaxRow));
disp('Chromosome: ');
disp(Top1);
disp('Time Elapsed:');
disp(TimeEl);