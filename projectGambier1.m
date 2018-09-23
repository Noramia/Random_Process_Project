%function projectGambier1(stoc_b,stoc_g,Infect,numsims)
%Input Variables
%---------------------------------------------
InfectA = 10;
InfectB = 0;


popA=200; %total population
popB=27000;
T=500; %maximum time
numsims=5;

%rem=1; %remediation option

%Parameters
gammaA=1/3;
gammaB=1/2;
betaA=1.6*gammaA;
betaB=1.3*gammaB;
% 
counter = [];
Cost = [];
Cases = [];
TimeTot = [];
MaxInf = [];

stoc_b = 0; %STOCHASTIC VARIANCE FOR BETA
stoc_g = 0; %AND gamma


%----------------------------------

%-----------------------------------

%----------------------------------------------
%Actual Loop-----------------------------------
%----------------------------------------------

for j=1:numsims
    counter = [counter,j];
    t = 0;%set initial time
    ts = 0;%start time vector
    istate2A=InfectA; %no of infected individuals at any time
    istate1A=popA-istate2A; %no of susceptible individuals
    istate2B=InfectB; %no of infected individuals at any time
    istate1B=popB-istate2B; %no of susceptible individuals
    ss1A = istate1A;%start infectious vector
    ss2A = istate2A;%start susceptible vector
    ss1B = istate1B;%start infectious vector
    ss2B = istate2B;%start susceptible vector
    
%For Loop J, whilst the disease exists
    while t<T %for t less than max time
        %Compute rates
        rIA=(betaA*istate2A*istate1A)/(popA-1)+(betaA*istate2B*istate1A)/(popA+popB-1);
        rIB=(betaB*istate2A*istate1B)/(popA+popB-1)+(betaB*istate2B*istate1B)/(popB-1);
        rRA=gammaA*istate2A;
        rRB=gammaB*istate2B;
        
        %For both situations possible
        if istate2A<popA
            if istate2B<popB
                sum1=rIA+rIB+rRA+rRB;
            t = t + exprnd(1/sum1);
            chose=rand;
            opt1=rIA/(sum1);
            opt2=opt1+rIB/sum1;
            opt3=opt2+rRA/sum1;
            if chose<opt1
                istate1A=istate1A-1;
                istate2A=istate2A+1;
            else if chose<opt2
                istate1B=istate1B-1;
                istate2B=istate2B+1;
                else if chose<opt3
                        istate2A=istate2A-1;
                    else
                        istate2B=istate2B-1;
                    end
                end
            end
            else
            sum2=rIA+rRA+rRB;
            t = t + exprnd(1/sum2);
            chose=rand;
            opt1=rIA/sum2;
            opt2=opt1+rRA/sum2;
            if chose<opt1
                istate1A=istate1A-1;
                istate2A=istate2A+1;
            else if chose<opt2
                  istate2A=istate2A-1;
                else
                  istate2B=istate2B-1;
                end
            end
            end
        else if istate2B<popB
                sum3=rIB+rRA+rRB;
                t = t+exprnd(1/sum3);
                chose=rand;
                opt1=rIB/sum3;
                opt2=opt1+rRA/sum3;
                if chose<opt1
                istate1B=istate1B-1;
                istate2B=istate2B+1;
                else if chose<opt2
                        istate2A=istate2A-1;
                    else
                        istate2B=istate2B-1;
                    end
                end
            else
                sum4=rRA+rRB;
                t=t+exprnd(1/sum4);
                if rand<rRA/sum4
                   istate2A=istate2A-1;
                else
                    istate2B=istate2B-1;
                end
            end
        end

        %Grow necessary vectors
        ts = [ts; t];
        ss1A = [ss1A; istate1A];
        ss1B = [ss1B; istate1B];
        ss2A = [ss2A; istate2A];
        ss2B = [ss2B; istate2B];
        ss = [ss1A,ss1B,ss2A,ss2B];
    if istate2A<1
        if istate2B<1
        break
        end
    end
    end

    
    %Cases =[Cases,pop-ss1(end)];
    %TimeTot = [TimeTot,ts(end-1)];
    %MaxInf = [MaxInf,max(ss2)];
    
    
    
end
plot(ts,ss2A)
hold
plot(ts,ss2B)
%cases= total number of infectios through the entire simulation
%avCase = mean(Cases);
%sdvCase = std(Cases); %std= standard deviation
%medCas=median(Cases)

% The time it takes for the infectios to die out.
%avTime = mean(TimeTot);
%sdvTime = std(TimeTot);
%medTime= median(TimeTot)

%The Max number of infections at any one time.
%avMax = mean(MaxInf);
%sdvMax = std(MaxInf);
%medMax= median(MaxInf)

%avCost = mean(Cost);
%medCost = median(Cost)
%sdvCost = std(Cost);
        
    
%plot(counter,Cost)
%histogram(Cost,100) % the distribution of the data itself
 
% look up plotting two plots at the same time!