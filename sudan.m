%% Non-Linear Timeseries Analysis Project 2023 - Testing on Rain Data from Sudan
% Aron Paulsson and Torbjörn Onshage
clear all;
close all;
addpath('functions', 'data');
clc

%% Load data
load('proj23.mat')

mod_t = ElGeneina.rain_org_t(1:300);
mod_y = ElGeneina.rain_org(1:300);

val_t = ElGeneina.rain_org_t(301:382);
val_y = ElGeneina.rain_org(301:382);

tst_t = ElGeneina.rain_org_t(383:end);
tst_y = ElGeneina.rain_org(383:end);

figure(); hold on;
plot(mod_t, mod_y)
plot(val_t, val_y)
plot(tst_t, tst_y)
legend('Model Data', 'Validation Data', 'Test Data')
title('Data Selection')
hold off

%% Initial Identification 

plot_residual(mod_y, 'El Geneina Rainfall', 25, true, 0.05)

plot_nl_residual(mod_y, 'El Geneina Rainfall', 25, 0.05, true, 0.05)

figure();
bcNormPlot(mod_y)
title('Box-Cox Normplot')
grid on;

%% Examination of Variance

[theta, knots] = smooth_spline(mod(mod_t, 1), mod_y, 0.00001);

mod_yhat = cox_deBoor(mod(mod_t, 1), knots{1}) * theta;
mod_resid = mod_y - mod_yhat;
MSE = 1/length(mod_yhat)*sum(mod_resid.^2)

ts = linspace(0, 1, 100);
yhat = cox_deBoor(ts, knots{1}) * theta;
%yhat = 4.8 .* ( cos(pi.*(ts -0.62) ).^2);

figure(); hold on; 
scatter(mod(mod_t, 1), mod_y)
plot(ts, yhat)
title('Yearly Trend of Data')
legend('Data', 'Smoothed Estimate')
xlabel('Time of Year as a Fraction')
ylabel('Rainfall [mm/month]')

N = 11;
edges = linspace(0, 1, N);

mod_y_var = zeros(1, N-1);

for i = 1:N-1
    idx = edges(i) <= mod(mod_t, 1) & mod(mod_t, 1) < edges(i+1);
    
    %mod_y_var(i) = sum(mod_resid( idx ).^2);
    mod_y_var(i) = var( mod_y( idx ) );
end

edges = [edges(1), repelem(edges(2:end-1), 2), edges(end)];

mod_y_var = repelem(mod_y_var, 2);

figure(); hold on;
plot(edges, mod_y_var, 'k-')
plot(linspace(0,1,100), 9000.*sin(pi.*linspace(0,1,100)).^2)
title('Variance Variation Over the Year')
hold off;

N = 30;
edges = linspace(min(mod_y), max(mod_y), N);

mod_y_var = zeros(1, N-1);

for i = 1:N-1
    idx = [false; edges(i) <= mod_y(1:end-1) & mod_y(1:end-1) < edges(i+1)];
    
    %mod_y_var(i) = sum(mod_resid( idx ).^2);
    mod_y_var(i) = var( mod_y( idx ) );
end

edges = [edges(1), repelem(edges(2:end-1), 2), edges(end)];

mod_y_var = repelem(mod_y_var, 2);

figure(); hold on;
title('Variance Variation Over Previous Value')
plot(edges, mod_y_var, 'k-') 
hold off;

%% Go For It

% initParameters = [ 1.35, 0.75, 6.3, 3];
% 
% func_derivative = @(t, X, THETA) [
%     THETA(2) .* ( THETA(1) - X(:,1) ), THETA(3) .* ones(length(X(:,1)), 1)
% ];

A = 6;
c = 5;
w = 6.25;

initParameters = [A, c, w];

initStates = [0, -1, 0.02];

func_derivative = @(t, X, param) [
    -param(3).*X(:,2) , param(3).*X(:,1) , param(2).*param(3).*X(:,1) .* exp( param(2).*( 1 + X(:,2) ) ) .* param(1).*exp(-2.*param(2));
];


%transitionFunction = @(X, THETA) X + 1/12 .* func_derivative(X, THETA);

transitionFunction = @(X, THETA) RK4step(func_derivative, THETA, 0, X, 1/12);

observationFunction = @(X, THETA) X(:,3);

numStates = 3;
numParameters = 3;

[estimatedStates, estimatedParameters] = AIF_estimate(log(1+mod_y), transitionFunction, observationFunction, numStates, numParameters, initStates, initParameters);


%%

yhat = zeros(length(mod_y), 1);

for i = 1:length(mod_y) 
    states = reshape(estimatedStates(i, end, :), 1, length(estimatedStates(end, end, :)));
    
    yhat(i,:) = observationFunction( states, estimatedParameters(end,:) );
end

yhat_val = zeros(length(val_y), 1);

states = reshape(estimatedStates(end, end, :), 1, length(estimatedStates(end, end, :)));

for i = 1:length(val_y)
    
    states = transitionFunction(states, initParameters);% estimatedParameters(end,:)); 
    
    yhat_val(i) = observationFunction(states, initParameters); % estimatedParameters(end,:));

end

figure(); hold on;
plot(mod_t, log(1+mod_y))
plot(mod_t, yhat)
plot(val_t, log(1+val_y))
%plot(val_t, yhat_val)
hold off;

%% State "Noise"

estStateTran = zeros(length(mod_y), 2);

estStateTran(1,:) = reshape(estimatedStates(1, end, :), 1, length(estimatedStates(1, end, :)));

for i = 2:length(mod_y)
    estStateTran(i,:) = transitionFunction( estStateTran(i-1,:), initParameters);%estimatedParameters(end,:));
end

figure(); hold on;
plot(mod_t, estimatedStates(:, end, 1), 'b')
plot(mod_t, estimatedStates(:, end, 2), 'r')
plot(mod_t, estStateTran(:,1), '--b')
plot(mod_t, estStateTran(:,2), '--r')
hold off;


%%

t = ElGeneina.rain_org_t;
y = ElGeneina.rain_org;

% figure();
% plot(t, y)

% logshift = 0.001;
% z = log(y+logshift);
% 
% figure();
% plot(t, z)

a1 = -4;
a2 = 180;
a3 = 50;
a4 = 0.55;

f = @(x, t) a1.*x - 2*a2*a3.*(mod(t, 1)-a4).*exp(- a3.*(mod(t, 1)-a4).^2);

yhat = zeros(length(y),1);
yhat(1) = y(1);
for i = 2:length(y)
    yhat(i) = y(i-1) + 1/12 * f(y(i-1), t(i-1));
    
    yhat(i) = max(0, yhat(i));
end

% figure(); hold on;
% plot(t, z)
% plot(t, zhat)
% hold off;
% yzhat = exp(zhat)-logshift;

figure(); hold on;
plot(t, y)
plot(t, yhat)
%plot(t, yzhat)
hold off;

y_residual = y - yhat;
%yz_residual = y - yzhat;

figure(); hold on;
plot(y_residual)
%plot(yz_residual)
hold off;

MSE = sum(y_residual.^2)

%yz_MSE = sum(yz_residual.^2)

%% New Model

t = ElGeneina.rain_org_t;
y = ElGeneina.rain_org;

func_derivative = @(t, X, param) [
    -(param(1).*param(4)./param(3)) .* X(2,:) - 0.1 .* X(1,:) .* ( (X(1,:)./param(1)).^2 + (X(2,:)./param(3)).^2 - 1 );
     (param(3).*param(4)./param(1)) .* X(1,:) - 0.1 .* X(2,:) .* ( (X(1,:)./param(1)).^2 + (X(2,:)./param(3)).^2 - 1 )
];

dt = 1/12;

func_transition = @(X, t) [ RK4step(func_derivative, X(3:6,:), t*dt, X(1:2,:), dt ); X(3:6,:) ] ...
                            + [1, 1, 1/(t^(2/3)).*ones(4,1)] .* randn(numParticles);

func_observation = @(y, X) normpdf( y, X(4,:) * (exp( X(5,:) + X(6,:) ) - 1), 5);

param = [3; 0.01; 5.2; 6.2];

X = zeros( 2, length(y) );
X(:,1) = [-0.6; -4.7];

sigma = 0.2;

for idx = 2:length(X)
    X(:,idx) = RK4step(func_derivative, param, t(idx), X(:,idx-1), dt) + sigma * sqrt(dt) * randn(2, 1);
end

figure(); hold on;
plot(t, y)
plot(t, param(2).*( exp( param(3) + X(2,:)) - 1) )
hold off


%% Full Series Prediction

yhat = zeros(length(y),1);
yhat(1) = y(1);
for i = 2:length(yhat)
    yhat(i) = yhat(i-1) + 1/12 * f(yhat(i-1), t(i-1));

    yhat(i) = max(0, yhat(i));
end

figure(); hold on;
plot(t, y)
plot(t, yhat)
hold off;

y_residual = y - yhat;
MSE = sum(y_residual.^2)

figure(); hold on;
plot(y_residual)
hold off;


%% Reconstruct Rain

upscale_factor = 2;

t_e = zeros( upscale_factor*(length(y)-1)+1,1 );

for idx_e = 1:length(t_e)-1
    idx = 1 + ( idx_e - 1 - mod(idx_e-1, upscale_factor) ) / upscale_factor;
    
    t_e(idx_e) = t(idx) + mod(idx_e-1, upscale_factor) / upscale_factor * (t(idx+1) - t(idx));
end
t_e(end) = t(end);

y_e = zeros( upscale_factor*(length(y)-1)+1,1 );

for idx = 1:length(y)-1
    idxs_e = ( 1+upscale_factor*(idx-1) ):( upscale_factor*idx );
    
    pred_left = zeros(upscale_factor,1); pred_left(1) = y(idx);
    pred_right = zeros(upscale_factor,1); pred_right(end) = y(idx+1);
    
    for j = 1:upscale_factor-1
        pred_left(1+j) = pred_left(j) + 1/(upscale_factor*12) .* f(pred_left(j), t_e(upscale_factor*(idx-1)+j));
        pred_right(end-j) = pred_left(end-j+1) + 1/(upscale_factor*12) .* f(pred_right(end-j+1), t_e(upscale_factor*idx-j+1));
    end
    
    pred_left = pred_left(2:end);
    pred_right = pred_right(1:end-1);
    
    y_e(idxs_e) = [y(idx); pred_left];% + pred_right) / 2];
end
y_e(end) = y(end);

y_e(y_e < 0) = 0;

figure(); hold on;
plot(t_e, y_e)
scatter(t, y)