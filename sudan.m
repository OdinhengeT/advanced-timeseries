%% Non-Linear Timeseries Analysis Project 2023 - Testing on Power Price Data
% Aron Paulsson and Torbjörn Onshage
clear all;
close all;
addpath('functions', 'data');
clc

%% Load data
load('proj23.mat')

%%

[theta, knots] = smooth_spline([100 .* mod(Kassala.rain_t(1:end-1), 1), Kassala.rain(1:end-1)], Kassala.rain(2:end), 0.01);

ts = linspace(0, 100, 100);
ys = linspace(0, 240, 100);

[X, Y] = meshgrid(ts, ys);

Yhat = 0.*X;

for i = 1:100
    Yhat(:,i) = [cox_deBoor(X(:,i), knots{1}), cox_deBoor(Y(:,i), knots{2})] * theta;
end

yhat = [cox_deBoor(100 .* mod(Kassala.rain_t(1:end-1), 1), knots{1}), cox_deBoor( Kassala.rain(1:end-1), knots{2})] * theta;
resid =  Kassala.rain(2:end) - yhat;
MSE = 1/length(yhat)*sum(resid.^2)

figure(); hold on;
scatter3(100 .* mod(Kassala.rain_t(1:end-1), 1), Kassala.rain(1:end-1), Kassala.rain(2:end))
surf(X, Y, Yhat)

%%

t = ElGeneina.rain_org_t;
y = log( ElGeneina.rain_org + 1 );

firsts_idx = find(conv(y > 0, [1, -1]) > 0); firsts_idx = firsts_idx(3:end-3);

histogram( mod( t(firsts_idx), 1 ), 0:1/14:1)

figure(); 
plot( mod(t(firsts_idx(1:end-1)),1) , mod(t(firsts_idx(2:end)),1) )

%%

[theta, knots] = smooth_spline(mod(t, 1), y, 0.00001);

yhat = cox_deBoor(mod(t, 1), knots{1}) * theta;
resid = y - yhat;
MSE = 1/length(yhat)*sum(resid.^2)

ts = linspace(0, 1, 100);
yhat = cox_deBoor(ts, knots{1}) * theta;
%yhat = 4.8 .* ( cos(pi.*(ts -0.62) ).^2);

figure(); hold on; 
scatter(mod(t, 1), y)
plot(ts, yhat)

figure(); plot(t, resid)

% Conditional Variance

N = 10;
edges = linspace(0, 1, N);

figure(); hold on;

y_var = zeros(1, N-1);

for i = 1:N-1
    idx = edges(i) <= mod(t, 1) & mod(t, 1) < edges(i+1);
    
    y_var(i) = sum(resid( idx ).^2);
end

edges = [edges(1), repelem(edges(2:end-1), 2), edges(end)];

y_var = repelem(y_var, 2);

plot(edges, y_var, 'k-') 



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