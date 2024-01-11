function [estimated_states, estimated_parameters] = EIF_estimate(data, ObservationDensity, TransitionDensity, numStates, numParams)

% Initialization
numParticles = 1000;
particles = randn(numParticles, numStates+numParams);                       % Initialize state particles randomly
T = length(data);
estimated_states = zeros(T,numStates);
estimated_parameters = zeros(T,numParams);
Theta_F = zeros(T, numParams); 

% SISR Algorithm
for t = 1:T
    % Mutation
    particles = TransitionDensity(particles, t, numParticles);

    % Weights
    weights = ObservationDensity(data(t), particles, numParticles);         % Measurement likelihood
    
    % Resampling
    numParticles = round(1000*exp(1-t) + 200 + 100*t^(0.6));
    [particles, weights] = resampling(particles, weights, numParticles);

    % State estimation
    weights = weights ./ sum(weights);                                      % Normalize weights
    estimated_states(t,:) = sum(weights .* particles(:,1));                 % Simplified since weights are already normalized

    % Parameter Estimation
    Theta_F(t,:) = sum(weights .* particles(:,numStates+1:end));          
    estimated_parameters(t,:) = Theta_F(t,:); 
end
end