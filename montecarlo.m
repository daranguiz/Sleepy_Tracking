%% Monte Carlo Simulation
% To determine the average number of time units it takes for an object
% following a random walk path to leave a 2D grid with uniform movement of
% up to distance 1.

clear all;
close all;

%% Parameters
tic
sims = 1000;
dist = 10;
NUM_STATES = (2*dist + 1)^2;

P = zeros(NUM_STATES);
for n = 1:NUM_STATES
    % First catch the four corners. BL,TR,TL,BR
    if(n==1)
        P(n,n) = 1/9;
        P(n,2) = 1/9;
        P(n,1+sqrt(NUM_STATES)) = 1/9;
        P(n,2+sqrt(NUM_STATES)) = 1/9;
    elseif(n==NUM_STATES)
        P(n,n) = 1/9;
        P(n,n-1) = 1/9;
        P(n,n-sqrt(NUM_STATES)) = 1/9;
        P(n,n-1-sqrt(NUM_STATES)) = 1/9;
    elseif(n==sqrt(NUM_STATES))
        P(n,n) = 1/9;
        P(n,n-1) = 1/9;
        P(n,n+sqrt(NUM_STATES)) = 1/9;
        P(n,n-1+sqrt(NUM_STATES)) = 1/9;
    elseif(n==NUM_STATES+1-sqrt(NUM_STATES))
        P(n,n) = 1/9;
        P(n,n+1) = 1/9;
        P(n,n-sqrt(NUM_STATES)) = 1/9;
        P(n,n+1-sqrt(NUM_STATES)) = 1/9;
    % Then catch the sides: L,B,R,T
    elseif(n<sqrt(NUM_STATES))
        P(n,n) = 1/9;
        P(n,n+1) = 1/9;
        P(n,n-1) = 1/9;
        P(n,n+sqrt(NUM_STATES)) = 1/9;
        P(n,n+1+sqrt(NUM_STATES)) = 1/9;
        P(n,n-1+sqrt(NUM_STATES)) = 1/9;
    elseif(mod(n,sqrt(NUM_STATES))==1)
        P(n,n) = 1/9;
        P(n,n+sqrt(NUM_STATES)) = 1/9;
        P(n,n-sqrt(NUM_STATES)) = 1/9;
        P(n,n+1) = 1/9;
        P(n,n+1+sqrt(NUM_STATES)) = 1/9;
        P(n,n+1-sqrt(NUM_STATES)) = 1/9;
    elseif(n>NUM_STATES-sqrt(NUM_STATES))
        P(n,n) = 1/9;
        P(n,n+1) = 1/9;
        P(n,n-1) = 1/9;
        P(n,n-sqrt(NUM_STATES)) = 1/9;
        P(n,n-1-sqrt(NUM_STATES)) = 1/9;
        P(n,n+1-sqrt(NUM_STATES)) = 1/9;
    elseif(mod(n,sqrt(NUM_STATES))==0)
        P(n,n) = 1/9;
        P(n,n-sqrt(NUM_STATES)) = 1/9;
        P(n,n+sqrt(NUM_STATES)) = 1/9;
        P(n,n-1) = 1/9;
        P(n,n-1-sqrt(NUM_STATES)) = 1/9;
        P(n,n-1+sqrt(NUM_STATES)) = 1/9;
    % Now do the interior.
    else
        P(n,n) = 1/9;
        P(n,n+1) = 1/9;
        P(n,n-1) = 1/9;
        P(n,n+sqrt(NUM_STATES)) = 1/9;
        P(n,n-sqrt(NUM_STATES)) = 1/9;
        P(n,n+1+sqrt(NUM_STATES)) = 1/9;
        P(n,n-1+sqrt(NUM_STATES)) = 1/9;
        P(n,n+1-sqrt(NUM_STATES)) = 1/9;
        P(n,n-1-sqrt(NUM_STATES)) = 1/9;
    end
end

%% Variable initialization
b_k = (dist*(2*dist + 1) + dist + 1)*ones(1,sims);
count = 0;
flag = 0;
count_arr = zeros(1,sims);

%% Simulations
while(~flag)

    % Increment time
    count = count + 1;

    % Generate a new object position IFF the object isn't done
    for i=1:sims
        if(b_k(i) ~= NUM_STATES+1)
            temp = [cumsum(P(b_k(i),:),2) 1];
            temp = find(rand < temp);
            b_k(i) = temp(1);
        end
    end

    count_arr(b_k~=NUM_STATES+1) = count_arr(b_k~=NUM_STATES+1) + ones(1,length(count_arr(b_k~=NUM_STATES+1)));
    
    % If we're done
    if(~any(b_k~=NUM_STATES+1))
        flag=1;
    end
end

%% Analysis
avg = sum(count_arr)/sims
toc