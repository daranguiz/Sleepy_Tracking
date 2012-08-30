% This is the nonlearning approach for values of T, uses the greedy
% algorithm.

function T = genT( P )
%GENT Generates the matrix of expected tracking cost increase should a
%particular sensor NOT turn on next time step.

global NUM_STATES;
global NUM_SENSORS;
global AWAKE_COST;
SIM_RUNS = 200;

% Now compute the various tracking costs we need
%rand('state',1029);
%randn('state',345);
T = zeros(NUM_STATES,NUM_SENSORS);
for state=1:NUM_STATES
    state
    prior = P(state,:);

    r = ones(1,NUM_SENSORS);
    % Perform SIM_RUNS Monte Carlo simulations to compute T; 200 default.
    flag = false;
    while (~flag)
        hypcosts = zeros(NUM_SENSORS,1);
        basecost = 0;
        for i=1:SIM_RUNS

            % Generate a realization of the state at the next time step
            temp = [cumsum(prior) 1];
            temp = find(rand < temp);
            b = temp(1);

            % Detect if we're done.
            if (b==NUM_STATES+1)
               continue;          
            end

            % Generate the observations
            s = genobs(b,(1:NUM_SENSORS));

            % Generate the weights for generating the posterior
            weights = genweights(s,(1:NUM_SENSORS));

            % Toggle the state of each sensor and record the tracking cost
            for sensor = 1:NUM_SENSORS
                
                % One at a time turn off a single sensor
                tempr=r;
                tempr(sensor) = ~tempr(sensor);

                % Start with p vector as the prior
                p = prior;
                
                % Looks at the sensor row of weights (because that's the
                % sensor we're concerned with) and multiplies this row of
                % weights with the prior.
                p = p .* prod(weights(tempr==0,:),1);
                
                % Renormalize p vector (to keep it a prob. dist.)
                p = p / sum(p);

                % Find the max element of the posterior dist.
                [~,bhat] = max(p);
                
                % Accumulate the hypothetical tracking costs for the
                % sensor.  Averages the costs as it goes by dividing by the
                % number of simulations.
                hypcosts(sensor) =  hypcosts(sensor) + 1/SIM_RUNS*(bhat ~= b);
            end

            % Compute the baseline tracking cost
            % r gets updated every time we get to the greedy assumption
            % code below.  The greedy code turns off a single sensor every
            % time, making the product take values other than one
            % eventually.
            p = prior;
            p = p .* prod(weights(r==0,:),1);
            p = p / sum(p);
            [~,bhat] = max(p);
            basecost = basecost + 1/SIM_RUNS*( bhat ~= b);
        end
%         % Code for All Asleep Assumption
%         flag = true;
         % Code for Greedy Assumption - The sensor that causes the largest
         % decrease in expected tracking cost is added to the awake set
         % until any further reduction due to a single sensor is less than
         % the awake cost.
         [maxval,maxpos] = max( basecost - hypcosts );
         if (maxval >= AWAKE_COST)
             r(maxpos) = 0;
         else
             flag = true;
         end
    end
    T(state,:) = abs(hypcosts' - basecost);
end

end