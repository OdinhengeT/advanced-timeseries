function [estimatedStates, estimatedParameters] = IF22_estimate(...
    time,                                                       ...
    data,                                                       ...
    systemFunctions,                                            ...
    numStates,                                                  ...
    numParams,                                              ...
    numIterations,                                              ...
    numParticles,                                               ...
    initialState,                                               ...
    initialParam,                                               ...
    rangeState,                                                 ...
    rangeParam)

if nargin < 6 ; numIterations = 25; end
    if nargin < 7 ; numParticles = 5000; end
    if nargin < 8 ; initialState = zeros(1, numStates); end
    if nargin < 9 ; initialParam = randn(1, numParam); end
    if nargin < 10; rangeState = repmat([-inf; inf], 1, numStates); end
    if nargin < 11; rangeParam = repmat([-inf; inf], 1, numParam); end

    if ~isvector(time); error("AIF_estumate: input 'time' must be vector"); end
    if ~iscolumn(time); time = time'; end
    
    if ~max( size(data) ~= length(time) ); error("AIF_estumate: input 'time' and 'data' must be same length"); end
    if size(data, 1) ~= length(time); data = data'; end

    if ~isvector(systemFunctions); error("AIF_estimate: input 'systemFunctions' must be vector"); end
    if length(systemFunctions) < 2; error("AIF_estimate: input 'systemFunctions' must be atleast length 2, containing the transition and observation functions"); end
    
    if all( size(rangeState) == [numStates, 2] ); rangeState = rangeState'; end
    if all( size(rangeState) ~= [2, numStates] ); error("AID_estimate: input 'rangeState' must be a matrix with size 2 by 'numState'"); end
    
    if all( size(rangeParam) == [numParams, 2] ); rangeParam = rangeParam'; end
    if all( size(rangeParam) ~= [2, numParams] ); error("AID_estimate: input 'rangeParam' must be a matrix with size 2 by 'numParam'"); end

M = numIterations;
N = length(data);
J = numParticles;

transitionFunction = systemFunctions{1};
observationFunction = systemFunctions{2};

estimatedStates = zeros(N,M,numStates);
estimatedParameters = zeros(M,numParams);

states_particles = initialState.*ones(J, numStates);
param_particles = initialParam.*ones(J, numStates);
Sigma = eye(numParams);
Sigma0 = Sigma;

h = @(theta_F, sigma_m) normrnd(theta_F, sigma_m*Sigma);
f_init = @(x_0, theta_F_0, sigma) normrnd(transitionFunction(0, x_0, theta_F_0, 0), sigma);
f_tran = @(x_F,theta_P, sigma) normrnd(transitionFunction(0, x_F, theta_P, 0), sigma);
f_obs_pdf = @(y,x_P,theta_P,sigma) normpdf(y,observationFunction(0, x_P, theta_P), sigma);

sigma_states = 1;
sigma_observations = 2;

theta_F = zeros(N,1);
V_P = zeros(N,1);

for m = 1:M
    m
    tau_m = 1/m;
    sigma_m = 1/(m^(1.5));
    param_particles_F_0 = h(param_particles, tau_m*Sigma);
    param_particles_F = param_particles_F_0;
    states_particles_F = f_init(states_particles, param_particles_F, 1);
    for n = 1:N
        param_particles_P = h(param_particles_F, sigma_m*Sigma);
        states_particles_P = f_tran(states_particles_F, param_particles_P, 1);
        param_particles_P(param_particles_P > 1) = 1;
        param_particles_P(param_particles_P < -1) = -1;
        weights = f_obs_pdf(data(n), states_particles_P, param_particles_P, 1);
        weights = weights / sum(weights);
        if sum(isnan(weights)) > 0
                fprintf('Warning: isnan(weights) != 0 \n')
        end
        if sum(weights) == 0
            fprintf('Warning: sum(weights) = 0 \n')
            weights = 1/J*ones(J,1);
        end

        % Evaluate Sigma
        if n == 1
            Sigma_prime = ((states_particles_P - states_particles_F).*weights)'*(states_particles_P - states_particles_F) ./ (tau_m + sigma_m);
        else
            Sigma_prime = Sigma_prime + ((states_particles_P - states_particles_F).*weights)'*(states_particles_P - states_particles_F) ./ sigma_m;
        end

        resampledIndices = randsample(1:J, J, true, weights);
        param_particles_F = param_particles_P(resampledIndices);
        states_particles_F = states_particles_P(resampledIndices);
        estimated_state = sum(weights .* states_particles_F);
        estimatedStates(n,m,:) = estimated_state;

        theta_F(n) = sum(weights.*param_particles_F);
        if n == 1
            theta_F_0 = mean(param_particles_F);
            V_P(n) = ((param_particles_P - theta_F_0).*weights)'*(param_particles_P - theta_F_0);
            sum_theta = V_P(n) \ (theta_F(n) - theta_F_0);
        else
            V_P(n) = ((param_particles_P - theta_F(n-1)).*weights)'*(param_particles_P - theta_F(n-1));
            sum_theta = sum_theta + V_P(n) \ (theta_F(n) - theta_F(n-1)); 
        end
    end
    param_particles = param_particles_F;
    estimatedParameters(m,:) = mean(param_particles);
    % if m == 1
    %    estimatedParameters(m,:) = initialParam + 1/m * sum_theta;
    % else
    %    estimatedParameters(m,:) = estimatedParameters(m-1,:) + 1/m * sum_theta;
    % end
    % param_particles = estimatedParameters(m,:)*ones(J, numStates);

    state_residual = estimatedStates(2:end,m,:) - f_tran(estimatedStates(1:end-1,m,:) , estimatedParameters(m,:), 0);
    sigma_states = std(state_residual);
    observation_residual = data - observationFunction(0, estimatedStates(1:end,m,:), estimatedParameters(m,:));
    sigma_observations = std(observation_residual);

    % Sigma estimate
    epsilon = 10^-2;
    Sigma_bis = Sigma_prime * norm(Sigma0) / norm(Sigma_prime);
    Sigma = (1 - epsilon)*Sigma_bis + epsilon*Sigma0;
end
end
