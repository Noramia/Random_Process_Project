%Random Processes Project: Mount Gambier simulation
%SW, SM, 11/10/2018

%For the following code: Suffix A denotes the value is school population
%                        Suffix B denotes the value is wider population
%                        State 1 is for susceptible
%                        State 2 is infectious
%Input Variables
%---------------------------------------------
InfectA = 10;%Initial infectious population
InfectB = 0;


popA=200; %total population
popB=25212;
T=500; %maximum time
numsims=5;%number of simulations. Run time is significant for >100 simulations

%----------------------------------------------

%Parameters
gammaA=1/3;
gammaB=1/2;
betaA=1.6*gammaA;
betaB=1.3*gammaB;

% Creating storage vectors
counter = [];
CasesA = [];
CasesB = [];
TimeTot = [];
MaxInfA = [];
MaxInfB = [];



%----------------------------------------------
%Simulation Loop-----------------------------------
%----------------------------------------------

for j=1:numsims
    counter = [counter,j];
    t = 0;%set initial time
    ts = 0;%start time vector
    istate2A=InfectA; %no of infected individuals school
    istate1A=popA-istate2A; %no of susceptible individuals school
    istate2B=InfectB; %no of infected individuals mt gam
    istate1B=popB-istate2B; %no of susceptible individuals mt gam
    ss1A = istate1A;%start sus vector school
    ss2A = istate2A;%start inf vector school
    ss1B = istate1B;%start sus vector mt gam
    ss2B = istate2B;%start inf vector mt gam
    
%For Loop J, whilst the disease exists
    while t<T %for t less than max time
        %Compute rates
        rIA=(betaA*istate2A*istate1A)/(popA-1)+(betaA*istate2B*istate1A)/(popA+popB-1); %Infection in school
        rIB=(betaB*istate2A*istate1B)/(popA+popB-1)+(betaB*istate2B*istate1B)/(popB-1); %infection in gambier
        rRA=gammaA*istate2A; %Recovery in school
        rRB=gammaB*istate2B; %Recovery in gambier
        
        %For both situations possible
        if istate2A<popA
            if istate2B<popB
                sum1=rIA+rIB+rRA+rRB;
            t = t + exprnd(1/sum1); %time step
            chose=rand;
            opt1=rIA/(sum1);%Creates a CDF of all probabilities
            opt2=opt1+rIB/sum1;
            opt3=opt2+rRA/sum1;
            if chose<opt1 %uses a random number to determine which event occurs
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
            %If there is infection possible in school, but not in wider
            %population
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
            %if infection is only possible in wider population
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
                %if no more infection is possible
                sum4=rRA+rRB;
                t=t+exprnd(1/sum4);
                if rand<rRA/sum4
                   istate2A=istate2A-1;
                else
                    istate2B=istate2B-1;
                end
            end
        end

        %Store values in each state at each time step
        ts = [ts; t];
        ss1A = [ss1A; istate1A];
        ss1B = [ss1B; istate1B];
        ss2A = [ss2A; istate2A];
        ss2B = [ss2B; istate2B];
        ss = [ss1A,ss1B,ss2A,ss2B];
        %suscept school, sus mt gamb, inf school, inf mt gambier
    if istate2A<1
        if istate2B<1
        break
        end
    end
    end

% store output values for each simulation    
CasesA =[CasesA,(popA-ss1A(end))];%Cases in school
CasesB =[CasesB,(popB-ss1B(end))];%Cases in mount gambier
TimeTot = [TimeTot,ts(end-1)];%total lifespan
MaxInfA = [MaxInfA,(max(ss2A))];%peak infectiion school
MaxInfB = [MaxInfB,(max(ss2B))];%peak infection mount gambier
     
end

%cases= total number of infectios through the entire simulation

Case = [mean(CasesA), std(CasesA), median(CasesA); mean(CasesB), std(CasesB), median(CasesB)]
Cases1=CasesA+CasesB
subplot(2,2,1), histogram((Cases1),100), xlabel('Total number of cases of infection'), ylabel('frequency')
% The time it takes for the infectios to die out.
TimeTot
Time = [mean(TimeTot), std(TimeTot), median(TimeTot)]
subplot(2,2,2), histogram(TimeTot,100),xlabel('Total length of outbreak (days)'), ylabel('frequency')
%The Max number of infections at any one time.
%MaxInfA
%MaxInfB
MaxInf = [mean(MaxInfA), std(MaxInfA), median(MaxInfA); mean(MaxInfB), std(MaxInfB), median(MaxInfB)]
MaxInfs= MaxInfA+MaxInfB
subplot(2,2,3), histogram((MaxInfs),100), xlabel('Max number of cases of infection'), ylabel('frequency')
