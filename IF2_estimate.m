function [estimatedStates, estimatedParameters] = IF2_estimate(data, transitionFunction, observationFunction, numStates, numParameters)
M = 100;
N = length(data);
J = 10000;

estimatedStates = zeros(N,M,numStates);
estimatedParameters = zeros(M,numParameters);

param_particles = zeros(J, numParameters);
states_particles = zeros(J, numStates);
Sigma = eye(numParameters);

h = @(theta_F, sigma_m) normrnd(theta_F, sigma_m*Sigma);
f_init = @(x_0, theta_F_0, sigma) normrnd(transitionFunction(x_0, theta_F_0), sigma);
f_tran = @(x_F,theta_P, sigma) normrnd(transitionFunction(x_F, theta_P), sigma);
f_obs_pdf = @(y,x_P,theta_P,sigma) normpdf(y,observationFunction(x_P, theta_P), sigma);

for m = 1:M
    m
    sigma_m = 1/(m);
    param_particles_F = h(param_particles, sigma_m);
    states_particles_F = f_init(states_particles, param_particles_F, 1);
    
    for n = 1:N
        param_particles_P = h(param_particles_F, sigma_m);
        states_particles_P = f_tran(states_particles_F, param_particles_P, 1);
        weights = f_obs_pdf(data(n), states_particles_P, param_particles_P, 1);
        weights = weights / sum(weights);
        resampledIndices = randsample(1:J, J, true, weights);
        param_particles_F = param_particles_P(resampledIndices);
        states_particles_F = states_particles_P(resampledIndices);
        estimatedStates(n,m,:) = sum(weights .* states_particles_F);
    end
    param_particles = param_particles_F;
    state_particles = states_particles_F;
    estimatedParameters(m,:) = mean(param_particles);
end
end
