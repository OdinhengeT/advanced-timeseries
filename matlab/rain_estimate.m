clear all
close all
clc
load('proj23.mat')
data = Kassala.rain_org;
%%
transitionFunction_cont = @(t0, states, params) [
    params(1)*(params(2) - states(1));
    params(3);
    exp(states(1)*(1+sin(states(2))) - states(3)) - 1/2];

observationFunction_cont = @(t0, states, params) params(4)*exp(states)-params(5)-1;

%%
clc
transitionFunction = @(states, params) RK4step(transitionFunction_cont, params, 0, states, 1/12);
observationFunction = @(states, params) RK4step(observationFunction_cont, params, 0, states, 1/12);

numStates = 3;
numParams = 5;

[states, params] = AIF_estimate(log(data + 1), transitionFunction, observationFunction, numStates, numParams);

%%
figure()
plot(exp(observationFunction_cont(0,states(:,end,3), params) - 1))
hold on
plot(data)
ylim([-10 300])
