%function EstProb = project1(pop, istate2, rem)

%Input Variables
%---------------------------------------------
pop=200; %total population
T=100; %maximum time

istate2=100; %no of infected individuals at any time
istate1=pop-istate2; %no of susceptible individuals

rem=0; %remediation option

%Parameters
gamma=1/3; 
beta=1.6*gamma;
mask=1;
%----------------------------------

%Remediation Alteration
%----------------------------------
if rem == 1
    gamma = 1/2;
end
if rem ==2
    mask=0.75;
end
if rem ==3;
    mask=0.75;
    gamma=1/2;
end
%-----------------------------------

    t = 0;%set initial time
    ts = 0;%start time vector
    ss2 = istate2;%start infectious vector
    ss1 = istate1;%start susceptible vector

    while t<T %for t less than max time
        rate1=(beta*istate2*istate1)/(pop-1);
        rate2=gamma*istate2;
        if istate2<200
            t = t + exprnd(1/(rate1+rate2));
            if rand<mask*(rate1/(rate1+rate2))
                istate1=istate1-1;
                istate2=istate2+1;
            else
                istate2=istate2-1;
            end
        else
            t = t + exprnd(rate2);
            istate2 = istate2-1;
        end
        ts = [ts; t];
        ss2 = [ss2; istate2];
        ss1 = [ss1; istate1];
    if istate2<1
        break
    end
    end

figure
plot(ts,ss1)
hold
plot(ts,ss2)
title('Plot of disease history')
xlabel('Time in days')
ylabel('Number of people')
legend('num of susceptible','num of infectious')