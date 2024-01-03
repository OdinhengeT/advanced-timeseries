function [estimatedStates, estimatedParameters] = AIF_estimate(data, transitionFunction, observationFunction, numStates, numParameters)
% Initialize
M = 25;
N = length(data);
J = 5000;

theta_ag = zeros(M,numParameters);
theta_ga = theta_ag;
theta = zeros(M,numParameters);

theta_n = zeros(N,numParameters);

theta_ag_0 = rand(1,numParameters);
theta_ga_0 = theta_ag_0;
theta_0 = zeros(1,numParameters);

a0 = 0.95^(0.02);
C = 1;

estimatedStates = zeros(N,M,numStates);
estimatedParameters = zeros(M,numParameters);
S = zeros(M,numParameters);

states_particles = randn(J, numStates);
Sigma = eye(numParameters);

L = 1/2;
beta = 5*1/(2*L);
lambda_0 = beta;

h = @(theta_F, sigma_m) mvnrnd(theta_F, sigma_m*Sigma, J);
f_tran = @(x_F,theta_P, sigma) transitionFunction(x_F, theta_P) + randn(1,numStates);     %normrnd(transitionFunction(x_F, theta_P), sigma);
f_obs_pdf = @(y,x_P,theta_P,sigma) normpdf(y,observationFunction(x_P(3), theta_P), sigma);

for m = 1:M
    fprintf('m = %i out of %i \n', m, M);
    sigma_m = 1/(m^(2/3));
    alpha = 2/(m+1);
    a = a0^(m-1);

    param_particles_F = h(theta_ga(m,:), (C*a*sigma_m));
    for j = 1:J
        states_particles_F(j,:) = f_tran(states_particles(j,:), param_particles_F(j,:), 1);                        % Initial states
    end
    
    if m == 1
        theta_ga(m,:) = (1-alpha)*theta_ag_0 + alpha*theta_0;
    else
        theta_ga(m,:) = (1-alpha)*theta_ag(m-1,:) + alpha*theta(m-1,:);
        lambda = (1+1/(m+1))*beta;
    end
    for n = 1:N
        if mod(n,N/100) == 0
            fprintf('%i percent of data filtered! \n', (n/N)*100);
        end
        param_particles_P = h(theta_ga(m,:), (a*sigma_m));                                           % Perturbe
        for j = 1:J
            states_particles_P(j,:) = f_tran(states_particles_F(j,:), param_particles_P(j,:), 1);                  % Simulate prediction particles
            weights(j,:) = f_obs_pdf(data(n), states_particles_P(j,:), param_particles_P(j,:), 1);                 % Evaluate weights
        end
        if sum(isnan(weights)) > 0
            fprintf('Warning: isnan(weights) != 0 \n')
        end
        if sum(weights) == 0
            fprintf('Warning: sum(weights) = 0 \n')
            weights = 1/J*ones(J,1);
        else
        weights = weights / sum(weights);                                                       % Normalize weights                   
        end              
        estimatedStates(n,m,:) = sum(weights .* states_particles_F);
        theta_n(n,:) = sum(weights .* param_particles_F);
        resampledIndices = randsample(1:J, J, true, weights);                                   % Resample
        param_particles_F = param_particles_P(resampledIndices,:);
        states_particles_F = states_particles_P(resampledIndices,:);
    end
    S(m,:) = (C^(-2)*a^(-2))*(Sigma\sum(theta_n - theta_ga(m,:))')'/(N+1);                       % Score function estimation
    if m == 1
        theta(m,:) = theta_0 + lambda_0*S(m,:);
        theta_ag(m,:) = theta_ga_0 + beta*S(m,:);
    else
        theta(m,:) = theta(m-1,:) + lambda*S(m,:);
        theta_ag(m,:) = theta_ga(m-1,:) + beta*S(m,:);
        lambda*S(m,:)
    end
    if abs(theta(m,:)) > 1
        theta(m,:) = 1*sign(theta(m,:));
    end
    estimatedParameters(m,:) = theta(m,:); 
end
end
