%% Non-Linear Timeseries Analysis Project 2023 - Testing on Sudan Data
% Aron Paulsson and Torbjörn Onshage
clear all;
close all;
addpath('functions', 'data');
clc

%% Load data
load('proj23.mat')

%% 3D Kassala rain

[theta, knots] = smooth_spline([100 .* mod(Kassala.rain_t(1:end-1), 1), Kassala.rain(1:end-1)], Kassala.rain(2:end), 0.01);

ts = linspace(0, 100, 100);
ys = linspace(0, 240, 100);

[X, Y] = meshgrid(ts, ys);

Yhat = 0.*X;

for i = 1:100
    Yhat(:,i) = [cox_deBoor(X(:,i), knots{1}), cox_deBoor(Y(:,i), knots{2})] * theta;
end

yhat = [cox_deBoor(100 .* mod(Kassala.rain_t(1:end-1), 1), knots{1}), cox_deBoor( Kassala.rain(1:end-1), knots{2})] * theta;
resid =  Kassala.rain(2:end) - yhat;
MSE = 1/length(yhat)*sum(resid.^2)

figure(); hold on;
scatter3(100 .* mod(Kassala.rain_t(1:end-1), 1), Kassala.rain(1:end-1), Kassala.rain(2:end))
surf(X, Y, Yhat)

%% Correlation between Monsune start between years?

figure(); hold on;
plot(t, y)
plot(t, 1+cos(8.5.*t))

firsts_idx = find(conv(y > 0, [1, -1]) > 0); firsts_idx = firsts_idx(3:end-3);

figure();
histogram( mod( t(firsts_idx), 1 ), 0:1/14:1)

figure(); 
plot( mod(t(firsts_idx(1:end-1)),1) , mod(t(firsts_idx(2:end)),1) )

%% INIT A MODEL
clear all
close all
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


%%

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

%% State "Noise"

figure(); hold on;
plot(mod_t, mod_state(:,1), 'b')
%plot(mod_t, mod_state(:,2), 'r')
plot(mod_t, estimatedStates(end, :, 1), '--b')
%plot(mod_t, estimatedStates(end, :, 2), '--r')
hold off;


