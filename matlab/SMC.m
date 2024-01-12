clear all
close all
clc
%% Generate Data %%
% Define parameters
T = 100;                                                                    % Number of time steps
numSims = 20;
% True parameter value
true_a = 0.7;

TransitionFunction_true = @(state) [true_a*state(:,1)] + [1] .* randn(1,1);
ObservationFunction_true = @(state) [state(:,1)] + 0.1*randn();

% Generate synthetic data
true_x = zeros(T,numSims);
Obs_Y = zeros(T,numSims);
for i =1:numSims
    obs_y = zeros(T,1);
    true_x(1,i) = randn(1,1);
    obs_y(1) = true_x(1,i) + randn;                                               % Initial observation
    
    for t = 2:T
        true_x(t,i) = TransitionFunction_true(true_x(t-1,i));                   % State transition equation
        obs_y(t) = ObservationFunction_true(true_x(t,i));                     % Observation equation
    end
    Obs_Y(:,i) = obs_y;
end

figure()
plot(Obs_Y)
grid on
title('Simulated realizations')
xlabel('iteration')
ylabel('value')

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
for i = 1:1
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

%%
data = Obs_Y(:,1);
transitionFunction = @(x, theta) theta(:,1).*(x(:,1));
observationFunction = @(x,theta) x(:,1);
numStates = 1;
numParams = 1;
figure()
for i = 1:1
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

% for i = 1:numSims
% figure()
% subplot(2,1,1)
% plot(a_res(25,i)*X_res(:,i))
% hold on
% plot(Obs_Y(:,i))
% legend('y pred','y')
% subplot(2,1,2)
% plot(Obs_Y(:,i) - a_res(25,i)*X_res(:,i))
% end
%%
clear a_res

data = Obs_Y(:,1);
transitionFunction = @(t, x, params, dt) [params(:,1).*x(:,1)];
observationFunction = @(t, x, params) x(:,1);
time = 1:length(Obs_Y(:,1));
systemFunctions = {transitionFunction, observationFunction};                                           
numIterations = 200;
numParticles = 1000;
initialState = 0;
initialParam = 0.7;
rangeState = [-inf;inf];
rangeParam = [-0.99; 0.99];
numStates = 1;
numParams = 1;

figure()
for i = 1:2
    [states, params] = IF22_estimate(time, ...
                                    Obs_Y(:,i), ...
                                    systemFunctions, ...
                                    numStates, ...
                                    numParams, ...
                                    numIterations, ...
                                    numParticles, ...
                                    initialState, ...
                                    initialParam, ...
                                    rangeState, ...
                                    rangeParam);

    X_res(:,i) = states(:,end);
    a_res(:,i) = params(:,1);

subplot(2,1,1)
plot(1:T, X_res(:,i))
hold on
plot(1:T, true_x(:,i))
title('SMC')
grid on

subplot(2,1,2)
plot(1:numIterations, a_res(:,i)-true_a*ones(numIterations,1))
hold on
plot(1:numIterations, true_a*ones(numIterations,1))
title('a residual')
grid on

MSE_x(i) = mean((X_res(:,i) - true_x(i)).^2);
MSE_a(i) = mean((a_res(:,i) - true_a).^2);
end
%%
figure()
boxplot(((a_res([1, 3, 7, 13, 21, 43, 73, 111, 157, 200],:)' - true_a)), [1, 3, 7, 13, 21, 43, 73, 111, 157, 200])
title('Convergence of a parameter over batches')
grid on
xlabel('batch')
ylabel('estimated a - true a')

% fprintf('IF2 MSE x\n')
% mean(MSE_x)
% fprintf('IF2 MSE a\n')
% mean(MSE_a)
% 
% 
for i = 1:numSims
figure()
subplot(2,1,1)
plot(a_res(end,i)*X_res(:,i))
hold on
plot(Obs_Y(:,i))
legend('y pred','y')
subplot(2,1,2)
plot(Obs_Y(:,i) - a_res(25,i)*X_res(:,i))
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

