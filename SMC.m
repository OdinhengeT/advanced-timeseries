clear all
close all
clc
%%
% Define parameters
T = 100; % Number of time steps
numParticles = 1000; % Number of particles

% True parameter value
true_a = 0.9;

% Generate synthetic data
true_x = zeros(T,1);
obs_y = zeros(T,1);
true_x(1) = randn; % Initial state
obs_y(1) = true_x(1) + randn; % Initial observation

for t = 2:T
    true_x(t) = true_a * true_x(t-1) + randn;           % State transition equation
    obs_y(t) = true_x(t) + randn;                    % Observation equation
end

% Initialize particles
particles_x = randn(numParticles, 1); % Initialize state particles randomly
particles_a = randn(numParticles, 1); % Initialize parameter particles randomly

% Initialize weights
weights = ones(numParticles, 1) / numParticles;

estimated_x = zeros(T,1);
estimated_a = zeros(T,1);

% SISR Algorithm
for t = 1:T
    % Prediction
    particles_x = particles_a .* particles_x + randn(numParticles, 1);
    
    % Update weights
    likelihood = normpdf(obs_y(t), particles_x, 1); % Measurement likelihood
    weights = weights .* likelihood;
    weights = weights / sum(weights);               % Normalize weights
    
    % Resampling
    resampledIndices = randsample(1:numParticles, numParticles, true, weights);
    particles_x = particles_x(resampledIndices);
    particles_a = particles_a(resampledIndices);
    weights = ones(numParticles, 1) / numParticles;
    
    % Parameter estimation (weighted mean)
    estimated_a(t) = sum(weights .* particles_a);
    
    % State estimation (weighted mean)
    estimated_x(t) = sum(weights .* particles_x);
end


%%
figure()
subplot(2,1,1)
plot(1:T,true_x)
hold on
plot(1:T,estimated_x)
legend('true x', 'estimated x')

subplot(2,1,2)
plot(1:T,true_a*ones(T,1))
hold on
plot(1:T,estimated_a)
legend('true a', 'estimated a')