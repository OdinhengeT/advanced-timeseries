function [estimated_states, estimated_parameters] = EIF_batch_estimate(data, observationFunction, transitionFunction, numStates, numParams)

f_tran = @(x_F,theta_P, sigma) transitionFunction(x_F, theta_P) + sigma.*randn(length(x_F),min(size(x_F)) + min(size(theta_P)));     %normrnd(transitionFunction(x_F, theta_P), sigma);
f_obs_pdf = @(y,x_P,theta_P,sigma) normpdf(y,observationFunction(x_P, theta_P), sigma);

% Initialization
J = 1000;
M = 10;
particles = randn(J, numStates+numParams);                       % Initialize state particles randomly
T = length(data);
Theta_F = zeros(T, numParams); 
estimated_states = zeros(T,numStates);
estimated_parameters = zeros(M,numParams);

Gamma = 1;
q = 1;
A = 0.05*T;

% SISR Algorithm
for m = 1:M
    for t = 1:T
        P = (1+q)*Gamma / (t + 1 + q);
        Q = P*q/(t+1+A);
        % Mutation
        particles = f_tran(particles(:,1:numStates), particles(:,numStates+1:end), Q);
    
        % Weights
        weights = f_obs_pdf(data(t), particles, J, Gamma);         % Measurement likelihood
        
        % Resampling
        J = round(1000*exp(1-t) + 200 + 100*t^(0.6));
        [particles, weights] = resampling(particles, weights, J);
    
        % State estimation
        weights = weights ./ sum(weights);                                      % Normalize weights
        estimated_states(t,:) = sum(weights .* particles(:,1:numStates));                 % Simplified since weights are already normalized
    
        % Parameter Estimation
        Theta_F(t,:) = sum(weights .* particles(:,numStates+1:end));           
    end
    estimated_parameters(m,:) = Theta_F(T,:);
end
end