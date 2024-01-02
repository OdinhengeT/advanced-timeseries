function [estimated_states, estimated_parameters] = EIF_estimate(data, ObservationDensity, TransitionDensity, numStates, numParams)
% Initialize particles
numParticles = 10000;
particles = randn(numParticles, numStates+numParams);                                            % Initialize state particles randomly
T = length(data);
estimated_states = zeros(T,numStates);
estimated_parameters = zeros(T,numParams);
Theta_F = zeros(T, numParams); 
V_P = zeros(T, numParams);
Phi = zeros(T, numParams); 
a = zeros(T, 1);
R = zeros(T, numParams); 
estimated_parameters( 1 ,:) = randn(1, numParams);
particles = TransitionDensity(particles, 1, numParticles);
weights = ObservationDensity(data(1), particles, numParticles);
weights = weights/sum(weights);
% SISR Algorithm
for t = 2:T
    % Mutation
    particles = TransitionDensity(particles, t, numParticles);
    V_P = (particles(:,numStates+1:end) - Theta_F(t-1,:))'*diag(weights)*(particles(:,numStates+1:end) - Theta_F(t-1,:));

    % Update weights
    weights = ObservationDensity(data(t), particles, numParticles);                       % Measurement likelihood
    
    % Selection and Resampling
    resampledIndices = randsample(1:numParticles, numParticles, true, weights);

    %numParticles = round(1000*exp(1-t) + 200 + 100*t^(0.6));
    if length(resampledIndices) > numParticles
        resampledIndices = resampledIndices(1:numParticles);
    elseif length(resampledIndices) < numParticles
        resampledIndices = randsample(resampledIndices, numParticles, true);
    end
    index = mod(resampledIndices,length(particles));
    index(index==0) = length(particles);
    particles = particles(index,:);
    weights = weights(index);

    % State estimation
    weights = weights ./ sum(weights);                                           % Normalize weights
    estimated_states(t,:) = sum(weights .* particles(:,1));                     % Simplified since weights are already normalized

    % Parameter Estimation
    Theta_F(t,:) = sum(weights .* particles(:,numStates+1:end));          
    Phi = t.*V_P; 
    %a = 1/(t + 1 + 0.05*T);
    a = 1/(t);
    R = V_P\(Theta_F(t,:) - Theta_F(t-1,:))';
    estimated_parameters(t,:) = estimated_parameters(t-1,:) + a*(Phi*R)'; 
    %estimated_parameters(t,:) = Theta_F(t,:);
end
end

