function [estimated_states] = SMC_estimate(data, ObservationDensity, TransitionDensity, numParticles)
if nargin < 4
    error("Error. Too few input arguments")
end

% Initialize particles
particles = randn(numParticles);                                            % Initialize state particles randomly
mean_a = 0;
T = length(data);
estimated_states = zeros(T,numParticles(2));

% SISR Algorithm
for t = 1:T
    % Mutation
    particles = TransitionDensity(particles, t);
    if t > 30                                                               % Use averaging after burn in pereiod
        mean_a = mean_a + (particles(:,2) - mean_a)/(t-30);
    end

    % Update weights
    weights = ObservationDensity(data(t), particles);                       % Measurement likelihood
    
    % Selection and Resampling
    resampledIndices = randsample(1:numParticles(1), numParticles(1), true, weights);
    particles = particles(resampledIndices,:);

    % State estimation
    weights = weights / sum(weights);                                       % Normalize weights
    estimated_states(t,:) = sum(weights .* particles);                      % Simplified since weights are already normalized
end
end

