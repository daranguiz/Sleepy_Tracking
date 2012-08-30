% In this 2D example the object can leave the network.
% It also assumes a Hamming tracking cost.

clear all;
close all;

%% Network type dependant initialization - 2D Network

global NUM_SENSORS;
NUM_SENSORS = 7^2;
global NUM_STATES;
NUM_STATES = (2*sqrt(NUM_SENSORS) - 1)^2;
global P;

% Create the transition probability matrix for a 2D network of
% states.  The target can move up to one space per time unit.
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

global sensx;
sensx = 2*floor(((1:NUM_SENSORS) - 1)/sqrt(NUM_SENSORS));
global statesx;
statesx = floor(((1:NUM_STATES) - 1)/sqrt(NUM_STATES));
global sensy;
sensy = 2*mod((1:NUM_SENSORS) - 1,sqrt(NUM_SENSORS));
global statesy;
statesy = mod((1:NUM_STATES) - 1,sqrt(NUM_STATES));
global positions;
positions = sensx;

% Start the target in the first position.
b_k = floor(NUM_STATES/2)+1;

%% Network type independant initialization
global AWAKE_COST;
AWAKE_COST = 0.025;
global T;
PAUSE_TIME = 0.1;   % 0.3 is pretty good for visual inspection.

% Start out p_k as knowing exactly where the target is.
%p_k = zeros(1,NUM_SENSORS+1);
%p_k(floor(NUM_SENSORS/2)+1) = 1;
% Start out p_k as a uniform distribution.
p_k = (1 / NUM_STATES ) * ones(1,NUM_STATES);

% Find the T vector of expected increase in tracking costs
T = genT(P);
%load('T7sensNetwork');
%load('T11sensNetwork');
%load('T13sensNetwork');

% Compute the initial sleeping policy
u = -1*ones(1,NUM_SENSORS);
r = zeros(1,NUM_SENSORS);
disp('DERP')
u(r==0) = sleeptimes( p_k, find(r==0), 1 );
disp('HERP')

%% Dynamic plot initialization
figure(1)
plot(sensx(:),sensy(:),'blacko');title('Sensor node map');
xlabel('X');ylabel('Y');

hold on;
plot(statesx(b_k),statesy(b_k),'ro');
pause(PAUSE_TIME);
hold off;

count = 0;
corr_count = 0;
time = 0;
ecost = 0;
tcost = 0;

%% Iterations
while (1)
    time = time + 1;
    r = r + u
    
    prevp = p_k;
    prevb = b_k;

    % Generate a new object position
    temp = [cumsum(P(b_k,:)) 1];
    temp = find(rand < temp);
    b_k = temp(1);

    % Terminate this simulation run if the object has left the network
    if (b_k == NUM_STATES+1)
        break;
    end;
    
    % Generate the observations
    s = -ones(1,NUM_SENSORS);
    s(r==0) = genobs(b_k,find(r==0));
    
    % Update the energy cost
    % New energy cost is the old one plus one for each of the awake sensors
    % in this time unit.
    ecost = ecost + sum(r==0);
    
    % Find average sleep times for the network
    avg_sleep_arr(time) = sum(r)/NUM_SENSORS;
    num_on_arr(time) = sum(r==0);
    
    % Compute a new posterior distribution
    prior = p_k*P;

    % Indexes of sensors that gave observations
    indices = find(r==0);
        
    weights = ones(NUM_SENSORS,NUM_STATES);
    
    % Posterior calculation
    weights(r==0,:) = genweights( s(indices), indices );
    p_k = prior .* prod(weights,1);
    p_k = p_k / sum(p_k);

    basecost = 1-max(p_k);

    % Generate some hypothetical observations for the sensors are not awake.
    % Since b is not known, we first generate b according to p, then
    % generate an observation, and update the weights
    temp = [cumsum(p_k) 1];
    temp = find(rand < temp);
    tempb = temp(1);

    s(r > 0) = genobs(tempb,find(r>0));
    weights(r>0,:) = genweights( s(r>0),find(r>0) );
    
    % Estimate the tracking cost associated with toggling the state of each sensor
    hypcosts = zeros(NUM_SENSORS,1);
    for sensor = 1:NUM_SENSORS
        tempr=r;
        if (tempr(sensor) == 0)
            tempr(sensor) = 1;
            
            % Recompute the posterior having thrown this observation away
            tempp = prior;
            tempp = tempp .* prod(weights(tempr==0,:),1);
            tempp = tempp / sum(tempp);

            hypcosts(sensor) =  1-max(tempp);
        else
            tempr(sensor) = 0;
            
            % Recompute the posterior having added the new observation
            tempp = p_k;
            tempp = tempp .* weights(sensor,:);
            tempp = tempp / sum(tempp);

            hypcosts(sensor) =  1-max(tempp);
        end
    end
    T = T - .01*prevp'*( prevp*T - (2*(r==0)-1).*(hypcosts'-basecost) );
    
    % Compute the estimate and update the tracking cost
    [temp,bhat] = max(p_k);
    if (bhat ~= b_k)
        tcost = tcost + 1;
    end
    tcostArr(time) = tcost;
    
    % Compute new sleep times
    u=-1*ones(1,NUM_SENSORS);
    % If any r (residual sleep times) are zero then let the u vector at
    % that sensor be updated with a new sleep time.
    if (any(r==0))
        u(r==0) = sleeptimes( p_k, find(r==0), 0 )
    end
    
    % Dynamic graph update
    figure(1)
    clf;
    plot(sensx(:),sensy(:),'blacko');title('Sensor node map');
    xlabel('X');ylabel('Y');
    hold on;
    plot(sensx(find(r==0)), sensy(find(r==0)), 'greeno');
    plot(statesx(b_k), statesy(b_k), 'redo');
    hold off;
    pause(PAUSE_TIME);
end;

%% Plot data
figure(2)
subplot(3,1,1);
plot(num_on_arr);xlabel('Time');ylabel('Number of sensors awake');
title('Number of Sensors Awake Vs. Time');
subplot(3,1,2);
plot(avg_sleep_arr);xlabel('Time');ylabel('Average sleep time per sensor');
title('Average Sleep Time');
subplot(3,1,3);
plot(tcostArr);xlabel('Time');ylabel('Tracking Error Counter');
title('Tracking Errors Vs. Time');
detection_rate = (time-tcost)/time
Avg_Energy_Use = sum(num_on_arr)/time
tracking_errors = tcost