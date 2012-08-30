function plotter(NUM_SENSORS)

% Plots the generated T values for a given state.  Assumes a grid 2D
% network.

close all;

% Loads some default T matrix if no parameters
if(nargin==0)
    NUM_SENSORS = 5^2;
end

% Check to see if it's something we can check
if(NUM_SENSORS == 5^2)
    genPlot(NUM_SENSORS);
elseif(NUM_SENSORS == 7^2)
    load('T7sensNetwork');
    NUM_STATES = (2*sqrt(NUM_SENSORS) - 1)^2;
    gridPlot(T,NUM_STATES,NUM_SENSORS);
elseif(NUM_SENSORS == 11^2)
    load('T11sensNetwork');
    NUM_STATES = (2*sqrt(NUM_SENSORS) - 1)^2;
    gridPlot(T,NUM_STATES,NUM_SENSORS);
elseif(NUM_SENSORS == 13^2)
    load('T13sensNetwork');
    NUM_STATES = (2*sqrt(NUM_SENSORS) - 1)^2;
    gridPlot(T,NUM_STATES,NUM_SENSORS);
else
    disp('Matrix doesnt exist for this number of sensors');
    return;
end

end

function gridPlot(T,NUM_STATES,NUM_SENSORS)

state = 21;

% Assume we're in state given.
statex = floor((state - 1)/sqrt(NUM_STATES));
statex = statex/2;
statey = mod(state - 1,sqrt(NUM_STATES));
statey = statey/2;
valueVect = T(state,:);

% Normalize valueVect... because I want to
valueVect = valueVect/sum(valueVect);

% Turn valueVect into a matrix.
Z = zeros(sqrt(NUM_SENSORS));
counter = 1;
for i=1:sqrt(NUM_SENSORS)
    for j=1:sqrt(NUM_SENSORS)
        Z(j,i) = valueVect(counter);
        counter = counter + 1;
    end
end

figure(1)
hold on;
surf((1:sqrt(NUM_SENSORS))-1,(1:sqrt(NUM_SENSORS))-1,Z);
plot3(statex,statey,max(max(Z)),'redo');
hold off;
grid on;
xlabel('X');ylabel('Y');zlabel('Value');

end

function genPlot(NUM_SENSORS)

state = 45;

if(NUM_SENSORS == 25)
    % Use
    % "save('T49state25sensRand','T','sensx','sensy','statesx','statesy')"
    % to save the file needed from an appropriate genT command.
    load('T100state25sensRand');
else
    return;
end

% Assume we're in state given.
statex = statesx(state);
statey = statesy(state);
valueVect = T(state,:);

% Normalize valueVect... because I want to
valueVect = valueVect/sum(valueVect);

figure(1)
hold on;
plot3(sensx(:),sensy(:),valueVect(:),'blueo');
plot3(statex,statey,max(max(valueVect)),'redo');
plot3(statesx(:),statesy(:),min(min(valueVect))*ones(length(statesx),1),'blackx');
hold off;
grid on;
xlabel('X');ylabel('Y');zlabel('Value');

end