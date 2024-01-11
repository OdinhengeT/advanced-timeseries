%% Non-Linear Timeseries Analysis Project 2023 - Implementation of the AIF Estimation Algorithm
% Aron Paulsson and Torbjörn Onshage
function [estimatedState, estimatedParam] = AIF_estimate(...
    time,                                                      ...
    data,                                                      ...
    systemFunctions,                                           ...
    numState,                                                  ...
    numParam,                                                  ...
    numIterations,                                             ...
    numParticles,                                              ...
    initialState,                                              ...
    initialParam,                                              ...
    rangeState,                                                ...
    rangeParam                                                 ...
) 
    if nargin < 6 ; numIterations = 25; end
    if nargin < 7 ; numParticles = 5000; end
    if nargin < 8 ; initialState = zeros(1, numState); end
    if nargin < 9 ; initialParam = randn(1, numParam); end
    if nargin < 10; rangeState = repmat([-inf; inf], 1, numState); end
    if nargin < 11; rangeParam = repmat([-inf; inf], 1, numParam); end

    if ~isvector(time); error("AIF_estumate: input 'time' must be vector"); end
    if ~iscolumn(time); time = time'; end
    
    if ~max( size(data) ~= length(time) ); error("AIF_estumate: input 'time' and 'data' must be same length"); end
    if size(data, 1) ~= length(time); data = data'; end

    if ~isvector(systemFunctions); error("AIF_estimate: input 'systemFunctions' must be vector"); end
    if length(systemFunctions) < 2; error("AIF_estimate: input 'systemFunctions' must be atleast length 2, containing the transition and observation functions"); end
    
    if all( size(rangeState) == [numState, 2] ); rangeState = rangeState'; end
    if all( size(rangeState) ~= [2, numState] ); error("AID_estimate: input 'rangeState' must be a matrix with size 2 by 'numState'"); end
    
    if all( size(rangeParam) == [numParam, 2] ); rangeParam = rangeParam'; end
    if all( size(rangeParam) ~= [2, numParam] ); error("AID_estimate: input 'rangeParam' must be a matrix with size 2 by 'numParam'"); end
    
    N = length(time);
    M = numIterations;
    J = numParticles;

    theta_ag = zeros(M,numParam);
    theta_ga = theta_ag;
    theta = zeros(M,numParam);

    theta_n = zeros(N,numParam);

    theta_ag_0 = initialParam;
    theta_ga_0 = theta_ag_0;
    theta_0 = initialParam;

    estimatedState = zeros(M, N, numState);
    estimatedParam = zeros(M, numParam);
    S = zeros(M,numParam);
    
    state_particles_F = zeros(J, numState);
    %states_particles = repmat(initialStates, J, 1) + 0.01 .* randn(J, numStates);
    
    % Parameter uncertainty proportional to rangeParam    
    %Sigma = 10 .* diag( min(rangeParam(2,:) - rangeParam(1,:), 100) );
    Sigma = eye(numParam);
    
    a0 = 0.95^(0.02);
    C = 1;
    
    L = 1/2;
    beta = 1/(2*L);
    
    %f_perturb = @(theta_F, sigma_m) mvnrnd(theta_F, sigma_m*Sigma, J);
    %f_transition = @(x_F, theta_P, sigma) transitionFunction(x_F, theta_P) + randn(1, numStates)*diag(sigma);    
    %f_observationP = @(y, x_P, theta_P, sigma) normpdf(y, f_observation(x_P, theta_P), sigma);
    
    f_transition = systemFunctions{1};
    f_observation = systemFunctions{2};
 
    if length(systemFunctions) < 3
        f_estimate_initial_state = @(x, time, data, param) ( f_observation(time(1), x, param) - data(1,:) ).^2 ...
            + ( f_observation(time(1), f_transition(time(1), x, param, time(2)-time(1)), param) - data(2,:) ).^2;
    else
        f_estimate_initial_state = systemFunctions{3};
    end   
    
    optimization_options = optimoptions('fminunc','Display','none');

    for m = 1:M
        fprintf('m = %i out of %i \n', m, M);
        sigma_m = 1/(m^(2/3));
        alpha_m = 2/(m+1);
        a_m = a0^(m-1);
        
        param_particles_F = mvnrnd( theta_ga(m,:), C*a_m*sigma_m*Sigma, J ); %f_perturb( theta_ga(m,:), (C*a_m*sigma_m) );
              
        for j = 1:J
            state_particles_F(j,:) = fminunc(@(x) f_estimate_initial_state(x, time, data, param_particles_F(j,:)), initialState, optimization_options);
        end
        %state_particles_F = f_transition(states_particles, param_particles_F, sqrt(time(2)-time(1)).*[0.1] ); % HARDCODED

        if m == 1
            theta_ga(m,:) = (1-alpha_m)*theta_ag_0 + alpha_m*theta_0;
            lambda = beta;
        else
            theta_ga(m,:) = (1-alpha_m)*theta_ag(m-1,:) + alpha_m*theta(m-1,:);
            lambda = (1+1/(m-1))*beta; 
        end
        
        for n = 1:N
            
            % Perturb
            %param_particles_P = mvnrnd( theta_ga(m,:), a_m*sigma_m*Sigma, J ); %f_perturb( theta_ga(m,:), (a_m*sigma_m) );

            uBnd_low = diag( normcdf(rangeParam(1,:), theta_ga(m,:), a_m*sigma_m*Sigma) );
            uBnd_high = diag( normcdf(rangeParam(2,:), theta_ga(m,:), a_m*sigma_m*Sigma) );
            
            u_param_particles_P = repmat(uBnd_low, J, 1) + rand(J, numParam) * (uBnd_high - uBnd_low);
            
            param_particles_P = zeros(J, numParam);
            
            for p = 1:numParam
                param_particles_P(:,p) = norminv(u_param_particles_P(:,p), theta_ga(m,:), a_m*sigma_m*Sigma);
            end          
            
            % Simulate prediction particles
            dt = time(max(2, n)) - time(max(1,n-1));
            state_particles_P = f_transition(time(max(1,n-1)), state_particles_F, param_particles_P, dt ) ...
                              + 2.*randn(size(initialState)); % ADD NOISE sqrt(dt)*0.01
            %% HARDCODED
            
            % Evaluate weights
            weights = normpdf(data(n), f_observation(time(max(1,n-1)), state_particles_P, param_particles_P), 1); % HARDCODED
            
            % weights = f_observeP(data(n), state_particles_P, param_particles_P, 1);  % HARDCODED
            
            if sum(isnan(weights)) > 0
                fprintf('Warning: isnan(weights) != 0 \n')
            end
            if sum(weights) == 0
                fprintf('Warning: sum(weights) = 0 \n')
                weights = 1/J*ones(J,1);
            else
                % Normalize weights      
                weights = weights / sum(weights);             
            end  
            
            estimatedState(m, n, :) = sum(weights .* state_particles_F);
            theta_n(n, :) = sum(weights .* param_particles_F);
            
            % Resample
            resampledIndices = randsample(1:J, J, true, weights);
            param_particles_F = param_particles_P(resampledIndices,:);
            state_particles_F = state_particles_P(resampledIndices,:);
        end
        
        % Score function estimation
        S(m,:) = (C^(-2)*a_m^(-2))*(Sigma\sum(theta_n - theta_ga(m,:))')'/(N+1);
        
        if m == 1
            theta(m,:) = theta_0 + lambda*S(m,:);
            theta_ag(m,:) = theta_ga_0 + beta*S(m,:);
        else
            theta(m,:) = theta(m-1,:) + lambda*S(m,:);
            theta_ag(m,:) = theta_ga(m-1,:) + beta*S(m,:);
            lambda*S(m,:)
        end
        estimatedParam(m,:) = theta(m,:); 
    end
end
