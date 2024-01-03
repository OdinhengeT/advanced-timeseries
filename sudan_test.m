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


