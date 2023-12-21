clear all
close all
clc
%% Generate Data %%
% Define parameters
T = 100;                                                                    % Number of time steps

% True parameter value
true_a = -0.5;

% Generate synthetic data
true_x = zeros(T,1);
obs_y = zeros(T,1);
true_x(1) = randn;                                                          % Initial state
obs_y(1) = true_x(1) + randn;                                               % Initial observation

for t = 2:T
    true_x(t) = true_a * exp(true_x(t-1)) + randn;                               % State transition equation
    obs_y(t) = true_x(t) + randn;                                           % Observation equation
end

plot(obs_y)

%%
% Define Observation and Transition densities
numParticles = [1000 2];
%ObservationDensity = @(y, state) state(:,1) + randn(numParticles(1), 1);      % Observation density,  p(y_k | x_k)  

ObservationDensity = @(y,state) normpdf(y,state(:,1),1); 
TransitionDensity = @(state, t) [state(:,2).*exp(state(:,1)), state(:,2)] + [1, 1/t] .* randn(numParticles);  % Transition density,   q(x_k+1 | x_k)

states = SMC_estimate(obs_y, ObservationDensity, TransitionDensity, numParticles);

%%
figure()
subplot(3,1,1)
plot(1:T,true_x)
hold on
plot(1:T,states(:,1))
legend('true x', 'estimated x')

subplot(3,1,2)
plot(1:T,true_a*ones(T,1))
hold on
plot(1:T,states(:,2))
legend('true a', 'estimated a')

subplot(3,1,3)
plot(1:T,abs(true_a*ones(T,1)-states(:,2)))
legend('a error')
%% Estimate the State Space %%
numParticles = 100;                                                        % Number of particles

% Define Observation and Transition densities
sigma_e = 1;
p = @(y,x) normpdf(y,x,sigma_e);                                            % Observation density,  p(y_k | x_k)
q = @(a,x) a .* exp(x) + randn(numParticles, 1);                                 % Transition density,   q(x_k+1 | x_k)

% Initialize particles
particles_x = randn(numParticles, 1);                                       % Initialize state particles randomly
particles_a = randn(numParticles, 1);                                       % Initialize parameter particles randomly

% Initialize weights
%weights = ones(numParticles, 1) / numParticles;

estimated_x = zeros(T,1);
estimated_a = zeros(T,1);

% SISR Algorithm
for t = 1:T
    % Mutation
    %particles_x = particles_a .* particles_x + randn(numParticles, 1);
    particles_x = q(particles_a, particles_x);
    particles_a = particles_a + 1/t * randn(numParticles,1);
    
    
    % Update weights
    weights = p(obs_y(t), particles_x);                                     % Measurement likelihood
    weights = weights / sum(weights);                                       % Normalize weights
    
    % Selection and Resampling
    resampledIndices = randsample(1:numParticles, numParticles, true, weights);
    particles_x = particles_x(resampledIndices);
    particles_a = particles_a(resampledIndices);

    % Reset weights
    %weights = ones(numParticles, 1) / numParticles;
    
    % Parameter estimation
    estimated_a(t) = sum(weights .* particles_a);                           % Simplified since weights are already normalized
    
    % State estimation
    estimated_x(t) = sum(weights .* particles_x);                           % Simplified since weights are already normalized
end


%%
figure()
subplot(3,1,1)
plot(1:T,true_x)
hold on
plot(1:T,estimated_x)
legend('true x', 'estimated x')

subplot(3,1,2)
plot(1:T,true_a*ones(T,1))
hold on
plot(1:T,estimated_a)
legend('true a', 'estimated a')

subplot(3,1,3)
plot(1:T,abs(true_a*ones(T,1)-estimated_a))
legend('a error')