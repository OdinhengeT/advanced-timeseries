%% Non-Linear Timeseries Analysis Project 2023 - Testing of Non-Parametric Methods smooth_spline and loess
% Aron Paulsson and Torbjörn Onshage
clear all;
close all;
addpath('functions');
clc

%% Test: Cox-de Boor

x = linspace(0, 1, 500);
knots = [0,0, 0, 0, linspace(0, 1, 11) 1, 1, 1, 1];

B = cox_deBoor(x, knots,4, 0);

figure(); hold on;
for idx = 1:length(B(1,:))
    plot(x, B(:, idx))
end

%% Test: LOESS and Smoothing Splines
N = 250;

x = linspace(-4, 12, N);
y = -0.02*x.^3 + 0.23.*x.^2 + randn(1,N);

%%

x_eval = linspace(-4.5, 12.5, 2*N);

lambda = NaN;%0.01; % NaN if algorithm should select
[theta, knots, lambda] = smooth_spline(x, y, lambda);
yhat_ss = cox_deBoor(x_eval, knots{1}, 4)*theta;

alpha = NaN;%0.15; % NaN if algorithm should select
[yhat_loess, alpha] = loess(x, y, x_eval, alpha);

figure(); hold on;
title('LOESS and Smoothing Splines Fit')
plot(x_eval, yhat_loess, '-')
plot(x_eval, yhat_ss, '-')
scatter(x, y, '.k');
xlabel('X')
ylabel('Y')
legend(['LOESS: ', num2str(alpha)], ['Smoothing Splines: ', num2str(lambda)], 'Data')

%% Test: Smoothing Splines, Multivariable

N = 2500;
Ns = 100;

xs = linspace(-0.1, 1.1, Ns);

x1 = randn(N, 1); 
x2 = randn(N, 1); 

%y = 4*x1.^2.*(1-x1) + 0.5.*(x2-0.5) +  0.01 .* randn(N, 1);
y = 0.2.*sin(2.*pi.*x1.*x2) + 0.01 .* randn(N,1);

[theta, knots] = smooth_spline([x1, x2] , y, 0.000002, 4, [-0.1, -0.1; 1.1, 1.1]);

yhat = [cox_deBoor(x1, knots{1}), cox_deBoor(x2, knots{2})] * theta;

resid = y - yhat;
MSE = 1/length(yhat)*sum(resid.^2)

[X, Y] = meshgrid(xs, xs);

Yhat = 0.*X;

for i = 1:Ns
    Yhat(:,i) = [cox_deBoor(X(:,i), knots{1}), cox_deBoor(Y(:,i), knots{2})] * theta;
end

figure(); hold on;
scatter3(x1, x2, y)
surf(X, Y, Yhat)
xlabel("x")
ylabel("y")
zlabel("z")

%% another test

N = 1000;
Ns = 100;

xs = linspace(0, 1, Ns);

y = zeros(N, 1);

y(1) = 0;
y(2) = 0.2;

r = 4;

for i = 3:N
    y(i) = r * y(i-1) * (1 - y(i-1)) + 0.005*randn(1); %  * exp(-0.05*abs(y(i-2)))
    if y(i) < 0 || y(i) > 1
        y(i) = r * y(i-1) * (1 - y(i-1));
    end
end

[theta, knots] = smooth_spline([y(1:end-2),y(2:end-1)] , y(3:end), 0.00002, 4);

yhat = [cox_deBoor(y(1:end-2), knots{1}), cox_deBoor(y(2:end-1), knots{2})]*theta;

resid = y(3:end)-yhat;

MSE = 1/length(yhat)*sum(resid.^2)

% figure(); hold on;
% plot(xs, yhat)
% scatter(y(1:end-1), y(2:end))
% hold off;

[X, Y] = meshgrid(xs, xs);

Yhat = 0.*X;

for i = 1:Ns
    Yhat(:,i) = [cox_deBoor(X(:,i), knots{1}), cox_deBoor(Y(:,i), knots{2})]*theta;
end

figure(); hold on;
scatter3(y(1:end-2), y(2:end-1), y(3:end))
surf(X, Y, Yhat)
xlabel("y(t-2)")
ylabel("y(t-1)")
zlabel("Y(t)")


