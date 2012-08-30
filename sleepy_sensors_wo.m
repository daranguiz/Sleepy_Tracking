% In this 2D grid simulation the object can leave the network.

clear all;
close all;

%% Parameter definitions
NUM_SENSORS = 11^2;
STATES = (2*sqrt(NUM_SENSORS)-1)^2 + 1;
AWAKE_COST = 0.002;
ITERATIONS = 50;
PAUSE_TIME = 0.01;   % 0.3 is pretty good for visual inspection.

% Create the transition probability matrix for a 2-D square sensor array.
% Example layout with 16 sensors:
% 4-----8-----12----16
% -------------------
% 3-----7-----11----15
% -------------------
% 2-----6-----10----14
% -------------------
% 1-----5-----9-----13
P = zeros(NUM_SENSORS+1);
P(NUM_SENSORS+1,NUM_SENSORS+1) = 1;
for n = 1:NUM_SENSORS
    % First catch the four corners. BL,TR,TL,BR
    if(n==1)
        P(n,n) = 0.2;
        P(n,2) = 0.2;
        P(n,1+sqrt(NUM_SENSORS)) = 0.2;
        P(n,2+sqrt(NUM_SENSORS)) = 0.2;
        P(n,NUM_SENSORS+1) = 0.2;
    elseif(n==NUM_SENSORS)
        P(n,n) = 0.2;
        P(n,n-1) = 0.2;
        P(n,n-sqrt(NUM_SENSORS)) = 0.2;
        P(n,n-1-sqrt(NUM_SENSORS)) = 0.2;
        P(n,NUM_SENSORS+1) = 0.2;
    elseif(n==sqrt(NUM_SENSORS))
        P(n,n) = 0.2;
        P(n,n-1) = 0.2;
        P(n,n+sqrt(NUM_SENSORS)) = 0.2;
        P(n,n-1+sqrt(NUM_SENSORS)) = 0.2;
        P(n,NUM_SENSORS+1) = 0.2;
    elseif(n==NUM_SENSORS+1-sqrt(NUM_SENSORS))
        P(n,n) = 0.2;
        P(n,n+1) = 0.2;
        P(n,n-sqrt(NUM_SENSORS)) = 0.2;
        P(n,n+1-sqrt(NUM_SENSORS)) = 0.2;
        P(n,NUM_SENSORS+1) = 0.2;
    % Then catch the sides: L,B,R,T
    elseif(n<sqrt(NUM_SENSORS))
        P(n,n) = 1/7;
        P(n,n+1) = 1/7;
        P(n,n-1) = 1/7;
        P(n,n+sqrt(NUM_SENSORS)) = 1/7;
        P(n,n+1+sqrt(NUM_SENSORS)) = 1/7;
        P(n,n-1+sqrt(NUM_SENSORS)) = 1/7;
        P(n,NUM_SENSORS+1) = 1/7;
    elseif(mod(n,sqrt(NUM_SENSORS))==1)
        P(n,n) = 1/7;
        P(n,n+sqrt(NUM_SENSORS)) = 1/7;
        P(n,n-sqrt(NUM_SENSORS)) = 1/7;
        P(n,n+1) = 1/7;
        P(n,n+1+sqrt(NUM_SENSORS)) = 1/7;
        P(n,n+1-sqrt(NUM_SENSORS)) = 1/7;
        P(n,NUM_SENSORS+1) = 1/7;
    elseif(n>NUM_SENSORS-sqrt(NUM_SENSORS))
        P(n,n) = 1/7;
        P(n,n+1) = 1/7;
        P(n,n-1) = 1/7;
        P(n,n-sqrt(NUM_SENSORS)) = 1/7;
        P(n,n-1-sqrt(NUM_SENSORS)) = 1/7;
        P(n,n+1-sqrt(NUM_SENSORS)) = 1/7;
        P(n,NUM_SENSORS+1) = 1/7;
    elseif(mod(n,sqrt(NUM_SENSORS))==0)
        P(n,n) = 1/7;
        P(n,n-sqrt(NUM_SENSORS)) = 1/7;
        P(n,n+sqrt(NUM_SENSORS)) = 1/7;
        P(n,n-1) = 1/7;
        P(n,n-1-sqrt(NUM_SENSORS)) = 1/7;
        P(n,n-1+sqrt(NUM_SENSORS)) = 1/7;
        P(n,NUM_SENSORS+1) = 1/7;
    % Now do the interior.
    else
        P(n,n) = 1/9;
        P(n,n+1) = 1/9;
        P(n,n-1) = 1/9;
        P(n,n+sqrt(NUM_SENSORS)) = 1/9;
        P(n,n-sqrt(NUM_SENSORS)) = 1/9;
        P(n,n+1+sqrt(NUM_SENSORS)) = 1/9;
        P(n,n-1+sqrt(NUM_SENSORS)) = 1/9;
        P(n,n+1-sqrt(NUM_SENSORS)) = 1/9;
        P(n,n-1-sqrt(NUM_SENSORS)) = 1/9;
    end
end

%% Sensor initialization
sens(NUM_SENSORS).x = 0;
for n = 1:NUM_SENSORS
    sens(n).x = floor((n-1)/sqrt(NUM_SENSORS));
    sens(n).y = mod((n-1),sqrt(NUM_SENSORS));
    sens(n).sleep_time = 0;
end
% Start the target in the first position.
b_k = floor(NUM_SENSORS/2);
% Start out p_k as knowing exactly where the target is.
p_k = zeros(NUM_SENSORS+1,1);
p_k(floor(NUM_SENSORS/2)) = 1;
p_k = p_k';
num_on_arr = zeros(ITERATIONS,1);
avg_sleep_arr = zeros(ITERATIONS,1);
tracking_err_arr = zeros(ITERATIONS,1);

%% Dynamic plot initialization
figure(1)
%grid on;
for n = 1:NUM_SENSORS
    hold on;
    plot(sens(n).x,sens(n).y,'blacko');title('Sensor node map');
    xlabel('X');ylabel('Y');
    hold off;
end
hold on;
plot(sens(b_k).x,sens(b_k).y,'ro');
pause(PAUSE_TIME);
hold off;
count = 0;
corr_count = 0;

%% Iterations
for l=1:ITERATIONS
    % Sets the last node to b_k for use by the line plotter function.
    last = b_k;
    
    % Generate the e_b_k vector using current b_k.
    e_b_k = zeros(1,NUM_SENSORS+1);
    for n = 1:(NUM_SENSORS+1)
        if(n == b_k)
            e_b_k(n) = 1;
        end
    end

    % Find the new b_k for the round.
    b_k = e_b_k*P;

    % Determine a new position from the probability distribution.
    temp = rand;
    for n = 1:(NUM_SENSORS+1)
        temp = temp - b_k(n);
        if(temp<0)
            b_k = n;
            break;
        end
    end
    
    % Break if we're done.
    if(b_k==NUM_SENSORS+1)
        break;
    end

    % Reset found and num_on to 0 for the new iteration.
    found = 0;
    num_on = 0;
    total_sleep_time = 0;

    % Read from the sensors.
    for n=1:NUM_SENSORS
        if(sens(n).sleep_time > 0)
            % If the sensor is asleep, let it keep sleeping...
            hold on;
            plot(sens(n).x,sens(n).y,'blacko');
            hold off;
            total_sleep_time = total_sleep_time + sens(n).sleep_time;
        else
            % If the sensor is awake, take a reading and decide how long to
            % sleep from there.
            if(b_k == n)
                % If the target is in sensor n's range...
                found = n;
            else
                % If the target is not in the range of sensor n...
                % We know it's NOT at sensor n.
                p_k(n) = 0;
                % Renormalize the probability matrix to reflect this.
                p_k = p_k./sum(p_k);
            end
            num_on = num_on + 1;
            hold on;
            plot(sens(n).x,sens(n).y,'go');
            hold off;
        end
    end
    
    % Build num_on array for post analysis.
    num_on_arr(l) = num_on;
    avg_sleep_time = total_sleep_time/NUM_SENSORS;
    avg_sleep_arr(l) = avg_sleep_time;

    % If we saw the target this iteration update the target position
    % probability distribution to reflect this. This results in p_k+1.
    if(found>0)
        p_k = zeros(1,NUM_SENSORS+1);
        p_k(found) = 1;
        count = count + 1;
    end
    
    % Find predicted b_k (ie highest p_k value)
    [~, I] = max(p_k);
    if(b_k == I)
        corr_count = corr_count + 1;
        if(l==1)
            tracking_err_arr(l) = 0;
        else
            tracking_err_arr(l) = tracking_err_arr(l-1);
        end
    else
        if(l==1)
            tracking_err_arr(l) = 1;
        else
            tracking_err_arr(l) = tracking_err_arr(l-1) + 1;
        end
    end
    
    % Find the new sleep time. (eq 20 from the paper)
    for n = 1:NUM_SENSORS
        if(sens(n).sleep_time>0)
            sens(n).sleep_time = sens(n).sleep_time - 1;
        else
            u = 0;
            temp = p_k*P;
            
            % Find the RHS of eq 20.
            sum_temp = 0;
            for i = 1:NUM_SENSORS
                if(temp(i) == 0)
                    % Do nothing.
                else
                    sum_temp = sum_temp + AWAKE_COST*temp(i);
                end
            end
            
            while((temp(n) < sum_temp) && (u<ITERATIONS))
                u = u + 1;
                temp = p_k*(P^(u+1));
                
                % Find the RHS of eq 20.
                sum_temp = 0;
                for i = 1:NUM_SENSORS
                    if(temp(i) <= 0)
                        % Do nothing.
                    else
                        sum_temp = sum_temp + AWAKE_COST*temp(i);
                    end
                end
            end
            sens(n).sleep_time = u;
        end
    end
    
    % Update p_k to show what we think p_k+1 is.
    p_k = p_k*P;
    
    % Update dynamic plot
    hold on;
    if(b_k<NUM_SENSORS+1)
        plot(sens(b_k).x,sens(b_k).y,'ro');
        line2(sens(last).x,sens(last).y,sens(b_k).x,sens(b_k).y);
    end
    hold off;
    pause(PAUSE_TIME)
end

%% Plot data
figure(2)
subplot(3,1,1);
plot(num_on_arr);xlabel('Time');ylabel('Number of sensors awake');
title('Number of Sensors Awake Vs. Time');
subplot(3,1,2);
plot(avg_sleep_arr);xlabel('Time');ylabel('Average sleep time per sensor');
title('Average Sleep Time');
subplot(3,1,3);
plot(tracking_err_arr);xlabel('Time');ylabel('Tracking Error Counter');
title('Tracking Errors Vs. Time');
detection_rate = count/l
Avg_Energy_Use = sum(num_on_arr)/l
tracking_error = 1 - corr_count/l