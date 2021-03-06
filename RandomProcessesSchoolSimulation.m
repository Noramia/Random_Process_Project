%Random Processes Project: School simulation
%SW,SM, 11/10/2018


%Input Variables
%---------------------------------------------
pop=200; %total population
T=200; %maximum time
numsims=50;%number of simulations
Infect = 10;%initial infected population

rem=0; %remediation option
%Remediation options are 0=No remediation
%                        1=Antiviral remediation
%                        2=Face Mask remediation
%                        3=Combined remediation

stoc_b=0.2;%laval of stochastic variance for infection
stoc_g=0.2;%level of stochastic variance for recovery
%---------------------------------------------


%Parameters (as provided in brief)
gamma=1/3; 
beta=1.6*gamma;
mask=1; %figure to potentially reduce infection rate with into of face masks (see if rem==2)
% 

%develop arrays to hold output results 
counter = [];
Cost = [];
Cases = [];
TimeTot = [];
MaxInf = [];

%----------------------------------

%Remediation Alteration
%----------------------------------
if rem == 1   % Just Antivirus
    gamma = 1/2;
end
if rem ==2 % Just the Mask
    mask=0.75;
end
if rem == 3  % both Mask and Antivirus
    gamma = 1/2;
    mask = 0.75;
end
%-----------------------------------

%----------------------------------------------
%Time Simulation-----------------------------------
%----------------------------------------------

for j=1:numsims
    counter = [counter,j];
    t = 0;%set initial time
    ts = 0;%start time vector
    istate2=Infect; %no of infected individuals at any time
    istate1=pop-istate2; %no of susceptible individuals
    ss2 = istate2;%start infectious vector
    ss1 = istate1;%start susceptible vector
    
%For Loop J, whilst the disease exists
    while t<T %for t less than max time
        %Compute rates
        rate1=(mask*(stoc_b*(2*rand-1)+beta)*istate2*istate1)/(pop-1);
        rate2=(stoc_g*(2*rand-1)+gamma)*istate2;
        %For both situations possible
        if istate2<200
            t = t + exprnd(1/(rate1+rate2));%Time to next event
            if rand<(rate1/(rate1+rate2))%probability of an infection
                istate1=istate1-1;%infection transition
                istate2=istate2+1;
            else
                istate2=istate2-1;%recovery transition
            end
        else
            %For only decrease possible
            t = t + exprnd(1/rate2);
            istate2 = istate2-1;%recovery transition
        end
        %Grow necessary vectors
        ts = [ts; t];
        ss2 = [ss2; istate2];
        ss1 = [ss1; istate1];
        ss = [ss1,ss2];
    if istate2<1
        break
    end
    end

    
    Cases =[Cases,pop-ss1(end)];%for given sim., store total cases
    TimeTot = [TimeTot,ts(end-1)];%for given sim, store lifespan
    MaxInf = [MaxInf,max(ss2)];%for given sim, store peak infection
    
    if rem == 1
    %for cost
    i=0;
    cost = 0;
    while i<T
        A = [ts(:,1)<(i+1) & ts(:,1)>=i];
        %times btwn d(i) and d(i+1)
        if sum(A)>0
          A2=zeros(size(A));
          for k=1:size(A)
          A2(k,k)=A(k);
          end
      B= A2*ss(:,2);
            %record infection no. ast each day i
            cost = cost + max(B)*2;
            % $2 per day for each infectious
        else cost = cost;
        end
        i=i+1;
    end
    Cost = [Cost, cost];
    end
    
    if rem == 2
    i=0;
   cost=0;
   while i<T
      A= [ts(:,1)<(i+1) & ts(:,1)>=i];
      if sum(A)>0
          %if any data points in this time add this day to cost
          cost=cost+ 0.5*pop;
      else cost = cost; 
      end
      i=i+1;
   end
   Cost= [Cost, cost];
    end

    if rem == 3
    %for cost
    i=0;
    cost = 0;
    while i<T
        A = [ts(:,1)<(i+1) & ts(:,1)>=i];
        if sum(A)>0
          A2=zeros(size(A));
          for k=1:size(A)
          A2(k,k)=A(k);
          end
      B= A2*ss(:,2);
            %record infection no. ast each day i
            cost = cost + max(B)*2 + pop*(0.5);
            % $2 per day for each infectious
        else cost = cost;
        end
        i=i+1;
    end
    Cost = [Cost, cost];
    end
    
        
        
end
%cases= total number of infectios through the entire simulation
avCase = mean(Cases); %Mean of total cases
sdvCase = std(Cases); %standard deviation of total cases
medCase=median(Cases); %median of total cases

%% The time it takes for the infectios to die out.
avTime = mean(TimeTot);
sdvTime = std(TimeTot);
medTime= median(TimeTot);

%The Max number of infections at any one time.
avMax = mean(MaxInf);
sdvMax = std(MaxInf);
medMax= median(MaxInf);

%Total cost
avCost = mean(Cost);
medCost = median(Cost);
sdvCost = std(Cost);
    
    
