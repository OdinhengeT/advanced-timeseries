function [estimated_states] = SMC_estimate_states(data, observationFunction, transitionFunction, numStates, params)
% Initialize particles
J = 1000;
particles = randn(J);                                            % Initialize state particles randomly
T = length(data);
estimated_states = zeros(T,numStates);

f_tran = @(x, sigma) normrnd(transitionFunction(x, -0.5), sigma);
f_obs_pdf = @(y,x,sigma) normpdf(y,observationFunction(x, -0.5), sigma);

% SISR Algorithm
for t = 1:T
    % Mutation
    particles = f_tran(particles,1);

    % Update weights
    weights = f_obs_pdf(data(t), particles, 1);                       % Measurement likelihood
    
    % Selection and Resampling
    resampledIndices = randsample(1:J, J, true, weights);
    particles = particles(resampledIndices,:);

    % State estimation
    weights = weights / sum(weights);                                       % Normalize weights
    estimated_states(t,:) = sum(weights .* particles);                      % Simplified since weights are already normalized
end
end

