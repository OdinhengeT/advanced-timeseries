function [estimated_states] = SMC_estimate(data, ObservationDensity, TransitionDensity, numParticles)
if nargin < 4
    error("Error. Too few input arguments")
end

% Initialize particles
particles = randn(numParticles);                                            % Initialize state particles randomly

T = length(data);
estimated_states = zeros(T,numParticles(2));

% SISR Algorithm
for t = 1:T
    % Mutation
    B = [1 1/t];
    particles = TransitionDensity(particles, B);
    
    % Update weights
    weights = ObservationDensity(data(t), particles);                       % Measurement likelihood
    
    % Selection and Resampling
    resampledIndices = randsample(1:numParticles(1), numParticles(1), true, weights);
    particles = particles(resampledIndices,:);

    % State estimation
    weights = weights / sum(weights);                                       % Normalize weights
    estimated_states(t,:) = sum(weights .* particles);                        % Simplified since weights are already normalized
end
end

