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

%% Go For It

A = 6;
c = 3;
w = 6.3;
mu = 0.3;

initParameters = [A, c, w, mu];
numParameters = length(initParameters);

initStates = [-1, -sqrt(3)]./2; %[-1, -1]./sqrt(2); %[0, -1];
numStates = length(initStates);

func_derivative = @(t, X, param) [
    -param(3).*X(:,2) - param(4).*X(:,1).*( X(:,1).^2 + X(:,2).^2 - 1 ), ...
     param(3).*X(:,1) - param(4).*X(:,2).*( X(:,1).^2 + X(:,2).^2 - 1 )
];

transitionFunction = @(X, THETA) RK4step(func_derivative, THETA, 0, X, 1/12);

observationFunction = @(X, THETA) THETA(1)*exp(-2.*THETA(2)) .* exp( THETA(2) .* (1+X(:,2)) );

mod_t = 0:1/12:100;
mod_y = zeros(length(mod_t), 1);
val_t = 100:1/12:150;
val_y = zeros(length(val_t), 1);

states = initStates;

for i = 1:length(mod_y)
    states = transitionFunction(states, initParameters);% + sqrt(1/12).*0.01.*randn(1,3); 
    mod_y(i) = observationFunction(states, initParameters);
end
for i = 1:length(val_y)
    states = transitionFunction(states, initParameters);% + sqrt(1/12).*0.01.*randn(1,3); 
    val_y(i) = observationFunction(states, initParameters);
end

[estimatedStates, estimatedParameters] = AIF_estimate(mod_t, mod_y, transitionFunction, observationFunction, numStates, numParameters, initStates, initParameters);


%%

yhat = zeros(length(mod_y), 1);

for i = 1:length(mod_y) 
    states = reshape(estimatedStates(i, end, :), 1, length(estimatedStates(end, end, :)));
    
    yhat(i,:) = observationFunction( states, estimatedParameters(end,:) );
end

yhat_val = zeros(length(val_y), 1);

states = reshape(estimatedStates(end, end, :), 1, length(estimatedStates(end, end, :)));

for i = 1:length(val_y)
    
    %states = transitionFunction(states, initParameters); 
    %yhat_val(i) = observationFunction(states, initParameters);
    
    states = transitionFunction(states, estimatedParameters(end,:)); 
    yhat_val(i) = observationFunction(states, estimatedParameters(end,:));

end

figure(); hold on;
plot(mod_t, log(1+mod_y), 'b')
plot(mod_t, yhat, 'r')
plot(val_t, log(1+val_y), 'b')
plot(val_t, yhat_val, 'r')
hold off;

%% State "Noise"

estStateTran = zeros(length(mod_y), numStates);

estStateTran(1,:) = reshape(estimatedStates(1, end, :), 1, length(estimatedStates(1, end, :)));

for i = 2:length(mod_y)
    %estStateTran(i,:) = transitionFunction( estStateTran(i-1,:), initParameters);
    estStateTran(i,:) = transitionFunction( estStateTran(i-1,:), estimatedParameters(end,:));
end

figure(); hold on;
plot(mod_t, estimatedStates(:, end, 1), 'b')
plot(mod_t, estimatedStates(:, end, 2), 'r')
plot(mod_t, estStateTran(:,1), '--b')
plot(mod_t, estStateTran(:,2), '--r')
hold off;


