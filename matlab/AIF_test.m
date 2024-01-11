%% Non-Linear Timeseries Analysis Project 2023 - Testing of the AIF Estimation Algorithm
% Aron Paulsson and Torbj�rn Onshage
% INIT AR(1)-MODEL
clear all
close all
addpath('functions', 'data');
clc

a = 0.7;

initParam = [a];
rangeParam = [0, 0.9];
numParam = length(initParam);

initState = [3];
rangeState = [-inf; inf];
numState = length(initState);

transitionFunction = @(t, X, data, dt) data(1) .* X(:,1);
observationFunction = @(t, X, data) X(:,1);
f_estimate_initial_state = @(x, time, data, param) ...
    ( observationFunction(time(1), transitionFunction(time(1), x, param, time(2)-time(1)), param) - data(1,:) ).^2;

systemFunctions = {transitionFunction, observationFunction, f_estimate_initial_state};


%% INIT CIRCLE MODEL
clear all
close all
clc

w = 1;
mu = 0.1;

initParam = [w];
rangeParam = [0; 1];
numParam = length(initParam);

initState = [-3 0]; %[-1, -1]./sqrt(2); %[0, -1];
rangeState  = [-inf, -inf; inf, inf];
numState = length(initState);

func_derivative = @(t, X, param) [
    -param(1).*X(:,2) - mu.*X(:,1).*( X(:,1).^2 + X(:,2).^2 - 1 ), ...
     param(1).*X(:,1) - mu.*X(:,2).*( X(:,1).^2 + X(:,2).^2 - 1 )
];

transitionFunction = @(t, X, THETA, dt) RK4step(func_derivative, THETA, t, X, dt);
observationFunction = @(t, X, THETA) X(:,2);
f_estimate_initial_state = @(x, time, data, param) ...
    ( observationFunction(time(1), transitionFunction(time(1), x, param, time(2)-time(1)), param) - data(1,:) ).^2;

systemFunctions = {transitionFunction, observationFunction, f_estimate_initial_state};

%% SIMULATE DATA

mod_t = (0:1/12:15)';
mod_state = zeros(length(mod_t), length(initState));
mod_y = zeros(length(mod_t), 1);

mod_state(1,:) = initState;
mod_y(1) = observationFunction(mod_t(1), mod_state(1,:), initParam);

for i = 2:length(mod_y)
    dt = mod_t(i) - mod_t(i-1);
    mod_state(i,:) = transitionFunction(mod_t(i), mod_state(i-1,:), initParam, dt) + 2.*randn(1,length(initState)); 
    mod_y(i) = observationFunction(mod_t(i), mod_state(i,:), initParam) + randn(1,1);
end

val_t = mod_t( ceil(length(mod_t)*2/3):end, : );
val_state = mod_state( ceil(length(mod_t)*2/3):end, : );
val_y = mod_y( ceil(length(mod_t)*2/3):end , : );

mod_t = mod_t( 1:floor(length(mod_t)*2/3), : );
mod_state = mod_state( 1:length(mod_t), : );
mod_y = mod_y( 1:length(mod_t), : );

%% ESTIMATE PARAM

[estimatedStates, estimatedParameters] = AIF_estimate(...
    mod_t,                                            ...
    mod_y,                                            ...
    systemFunctions,                                  ...
    numState,                                         ...
    numParam,                                         ...
    50,                                               ...
    7500,                                             ...
    initState,                                        ...
    initParam,                                        ...
    rangeState,                                       ...
    rangeParam                                        ...
);


%% Compare Final Predictions to Data

mod_yhat = zeros(length(mod_y), 1);

for i = 1:length(mod_y) 
    state = reshape(estimatedStates(end, i, :), 1, length(estimatedStates(end, end, :)));
    
    mod_yhat(i,:) = observationFunction(mod_t(i), state, estimatedParameters(end,:) );
end

val_yhat = zeros(length(val_y), 1);

state = reshape(estimatedStates(end, end, :), 1, length(estimatedStates(end, end, :)));

for i = 1:length(val_y)
    
    %states = transitionFunction(states, initParameters); 
    %yhat_val(i) = observationFunction(states, initParameters);
    
    dt = val_t( max(2,i) ) - val_t( max(1,i-1));
    state = transitionFunction(val_t(i), state, estimatedParameters(end,:), dt); 
    val_yhat(i) = observationFunction(val_t(i),state, estimatedParameters(end,:));

end

figure(); hold on;
plot(mod_t, mod_y, 'b')
plot(mod_t, mod_yhat, 'r')
plot(val_t, val_y, 'b')
plot(val_t, val_yhat, 'r')
hold off;

%% COMPARE ESTIMATED STATES TO TRUE STATES

figure(); hold on;
plot(mod_t, mod_state(:,1), 'b')
%plot(mod_t, mod_state(:,2), 'r')
plot(mod_t, estimatedStates(end, :, 1), '--b')
%plot(mod_t, estimatedStates(end, :, 2), '--r')
hold off;