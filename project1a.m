%function EstProb = project1a(pop, istate2, rem)
clear
clc
%Input Variables
%---------------------------------------------
pop=200; %total population
T=50; %maximum time
numsims=5000;

rem=1; %remediation option

%Parameters
gamma=1/3; 
beta=1.6*gamma;
mask=1;
% 
counter = [];
Cost = [];
Cases = [];
TimeTot = [];
MaxInf = [];

stoc_b = 0; %STOCHASTIC VARIANCE FOR BETA
stoc_g = 0; %AND gamma


%----------------------------------

%Remediation Alteration
%----------------------------------
if rem == 1
    gamma = 1/2;
end
if rem ==2
    mask=0.75;
end
if rem == 3
    gamma = 1/2;
    mask = 0.75;
end
%-----------------------------------

%----------------------------------------------
%Actual Loop-----------------------------------
%----------------------------------------------

for j=1:numsims
    counter = [counter,j];
    t = 0;%set initial time
    ts = 0;%start time vector
    istate2=10; %no of infected individuals at any time
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
            t = t + exprnd(1/(rate1+rate2));
            if rand<(rate1/(rate1+rate2))
                istate1=istate1-1;
                istate2=istate2+1;
            else
                istate2=istate2-1;
            end
        else
            %For only decrease possible
            t = t + exprnd(1/rate2);
            istate2 = istate2-1;
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

    
    Cases =[Cases,pop-ss1(end)];
    TimeTot = [TimeTot,ts(end-1)];
    MaxInf = [MaxInf,max(ss2)];
    
    if rem == 1
    %for cost
    i=0;
    cost = 0;
    X = [];
    while i<T
        A = [ts(:,1)<(i+1) & ts(:,1)>=i];
        %times btwn d(i) and d(i+1)
        if sum(A)>0
          A2=zeros(size(A));
          for k=1:size(A)
          A2(k,k)=A(k);
          end
      B= A2*ss(:,2);
      X= [X,max(B)];
            %record infection no. ast each day i
            cost = cost + max(B)*2;
            % $2 per day for each infectious
        else cost = cost;
        end
        i=i+1;
    end
    Cost = [Cost, cost];
    X=X;
    end
    
    if rem == 2
        
           %Estimating cost
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
        
        
end
    
avCase = mean(Cases);
sdvCase = std(Cases);

avTime = mean(TimeTot);
sdvTime = std(TimeTot);

avMax = mean(MaxInf);
sdvMax = std(MaxInf);

avCost = mean(Cost)
medCost = median(Cost)
sdvCost = std(Cost)
        
    
%plot(counter,Cost)
%histogram(Cost,100)


%Cases
%TimeTot
%MaxInf




%plot(ts,ss1)
%hold
%plot(ts,ss2)