%% Non-Linear Timeseries Analysis Project 2023 - Testing of the SMC Estimation Algorithm
% Aron Paulsson and Torbjörn Onshage
clear all
close all
clc
%% Generate Data %%
% Define parameters
T = 100;                                                                    % Number of time steps
numSims = 1;
% True parameter value
true_a = 0.2;

TransitionFunction_true = @(state) [true_a*state(:,1)] + [1] .* randn(1,1);
ObservationFunction_true = @(state) [state(:,1)] + 0.1*randn();

% Generate synthetic data
for i =1:numSims
    true_x = zeros(T,1);
    obs_y = zeros(T,1);
    Obs_Y = zeros(T,numSims);
    true_x(1,:) = randn(1,1);
    obs_y(1) = true_x(1) + randn;                                               % Initial observation
    
    for t = 3:T
        true_x(t,:) = TransitionFunction_true(true_x(t-1,:));                   % State transition equation
        obs_y(t) = ObservationFunction_true(true_x(t,:));                     % Observation equation
    end
    Obs_Y(:,i) = obs_y;
end

plot(obs_y)
hold on
plot(true_x(:,1))

%%
tic
% Define Observation and Transition densities
numParticles = [1000 3];
ObservationDensity_pdf = @(y,state) normpdf(y,state(:,1),0.1); 
TransitionDensity = @(state, t) [(state(:,2)+state(:,3)).*state(:,1), state(:,2), state(:,3)] + [1, 1, 1/t] .* randn(numParticles);  % Transition density,   q(x_k+1 | x_k)

X_res = zeros(T,numSims);
a_res = zeros(T,numSims);
figure()
for i = 1:numSims
    states = SMC_estimate(Obs_Y(:,i), ObservationDensity_pdf, TransitionDensity, numParticles);

    X_res(:,i) = states(:,1);
    a_res(:,i) = states(:,2);

subplot(2,1,1)
hold on
plot(1:T,X_res(:,i)-true_x(:,i))
plot(1:T,zeros(T,1))
legend('x error')
title('SMC')
subplot(2,1,2)
hold on
plot(1:T,a_res(:,i)-true_a(i))
plot(1:T,zeros(T,1))
legend('a error')
title('SMC')

MSE_x(i) = mean((X_res(:,i) - true_x(:,i)).^2);
MSE_a(i) = mean((a_res(:,i) - true_a(i)).^2);
end
%figure()
%boxplot(a_res(1:10:end,:)' - true_a, 1:10:100)
%title('SMC')

fprintf('SMC MSE X\n')
mean(MSE_x)
fprintf('SMC MSE a\n')
mean(MSE_a)
toc

for i = 1:numSims
figure()
subplot(4,1,1)
plot(a_res(end,i).*X_res(:,i))
hold on
plot(Obs_Y(:,i))
title('Y vs prediction')
legend('y pred','y')
subplot(4,1,2)
plot(Obs_Y(:,i) - a_res(end,i).*X_res(:,i))
title('Y resid')
subplot(4,1,3)
plot(true_a(i)*true_x(:,i))
hold on
plot(Obs_Y(:,i))
legend('y pred','y')
title('Y vs best possible')
subplot(4,1,4)
plot(Obs_Y(:,i) - true_a(i)*true_x(:,i))
title('best possible')
end
%%
tic
% Define Observation and Transition densities
numParticles = 1000;
numStates = 1;
numParams = 1;
ObservationDensity = @(y,state, numParticles) normpdf(y,state(:,1),1); 
TransitionDensity = @(state, t, numParticles) [state(:,2).*(state(:,1)), state(:,2)] + [1, 1/t] .* randn(numParticles, numStates+numParams);  % Transition density,   q(x_k+1 | x_k)
X_res = zeros(T,numSims);
a_res = zeros(T,numSims);
%b_res = zeros(T,numSims);

figure()
for i = 1:numSims
    [states, params] = EIF_estimate(Obs_Y(:,i), ObservationDensity, TransitionDensity, numStates, numParams);

    X_res(:,i) = states(:,1);
    a_res(:,i) = params(:,1);

subplot(2,1,1)
plot(1:T,states(:,1)-true_x(i))
hold on
plot(1:T,zeros(T,1))
legend('x error')
title('SMC')
subplot(2,1,2)
plot(1:T,params(:,1)-true_a(i))
hold on
plot(1:T,zeros(T,1))
legend('a error')
title('SMC')

MSE_x(i) = mean((X_res(:,i) - true_x(i)).^2);
MSE_a(i) = mean((a_res(:,i) - true_a(i)).^2);
end
figure()
boxplot(a_res(1:10:end,:)' - true_a, 1:10:100)
title('EIF')
fprintf('EIF MSE x\n')
mean(MSE_x)
fprintf('EIF MSE a\n')
mean(MSE_a)
toc
%%
tic
data = Obs_Y(:,1);
transitionFunction = @(x, theta) theta(:,1).*(x(:,1));
observationFunction = @(x,theta) x(:,1);
numStates = 1;
numParams = 1;
figure()
for i = 1:numSims
    i
    [states, params] = IF2_estimate(Obs_Y(:,i), transitionFunction, observationFunction, numStates, numParams);
    
    X_res(:,i) = states(:,end,1);
    a_res(:,i) = params(:,1);

subplot(2,1,1)
plot(1:T,states(:,1)-true_x(i))
hold on
plot(1:T,zeros(T,1))
legend('x error')
title('SMC')
subplot(2,1,2)
plot(1:T,params(:,1)-true_a(i))
hold on
plot(1:T,zeros(T,1))
legend('a error')
title('SMC')

MSE_x(i) = mean((X_res(:,i) - true_x(i)).^2);
MSE_a(i) = mean((a_res(:,i) - true_a(i)).^2);
end
figure()
boxplot(a_res(1:10:end,:)' - true_a, 1:10:100)
title('IF2')

fprintf('IF2 MSE x\n')
mean(MSE_x)
fprintf('IF2 MSE a\n')
mean(MSE_a)
toc

for i = 1:numSims
figure()
subplot(2,1,1)
plot(a_res(25,i)*X_res(:,i))
hold on
plot(Obs_Y(:,i))
legend('y pred','y')
subplot(2,1,2)
plot(Obs_Y(:,i) - a_res(25,i)*X_res(:,i))
end
%%
clear a_res
tic
data = Obs_Y(:,1);
transitionFunction = @(x, theta) [theta(:,1).*(x(:,1) + x(:,2)), x(:,2)];
observationFunction = @(x,theta) x(:,1);
numStates = 2;
numParams = 1;
figure()
for i = 1:numSims
    i
    [states, params] = AIF_estimate(Obs_Y(:,i), transitionFunction, observationFunction, numStates, numParams);
    
    X_res(:,i) = states(:,end,1);
    X2_res(:,i) = states(:,end,2);
    a_res(:,i) = params(:,1);

subplot(3,1,1)
plot(1:T,states(:,1))
hold on
plot(1:T,true_x(:,1))
legend('x_1')
title('AIF')
subplot(3,1,2)
plot(1:T,states(:,2))
hold on
plot(1:T,true_x(:,2))
legend('x_2')
title('AIF')
subplot(3,1,3)
plot(1:25,params(:,1))
hold on
plot(1:25,true_a(i)*ones(25,1))
legend('a')
title('AIF')
end
toc

for i = 1:numSims
figure()
subplot(2,1,1)
plot(a_res(25,i).*X_res(:,i).*X2_res(:,i) + X2_res(:,i))
hold on
plot(Obs_Y(:,i))
legend('y pred','y')
subplot(2,1,2)
plot(Obs_Y(:,i) - a_res(25,i).*X_res(:,i).*X2_res(:,i) + X2_res(:,i))
legend('y pred resid')
end

%%
transitionFunction = @(x, theta) [theta(:,1).*(x(:,1)), theta(:,1)];
observationFunction = @(x,theta) x(:,1);
numStates = 1;
numParams = 1;

[states, params] = EIF_batch_estimate(Obs_Y, observationFunction, transitionFunction, numStates, numParams);

%%
figure()
plot(states(:,1))
hold on
plot(true_x(:,1))

figure()
plot(true_a.*ones(10,1));
hold on
plot(params)

