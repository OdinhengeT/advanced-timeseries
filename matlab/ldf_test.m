%% Non-Linear Timeseries Analysis Project 2023 - Testing Example of LDF and PLDF
% Aron Paulsson and Torbjörn Onshage
clear all;
close all;
addpath('functions');
clc

%% Chaotic Logistic Map Test

N = 100;

y = zeros(N, 1);

y(1) = 0.29;

r = 4;

for i = 2:N
    y(i) = r * y(i-1) * (1 - y(i-1)); % + 0.05*randn(1);
    if y(i) < 0 || y(i) > 1
        y(i) = r * y(i-1) * (1 - y(i-1));
    end
end

lambda = NaN; %0.0001;

figure();
subplot(2, 3, 1); hold on;
    plot(y)
    title('Realization')
    xlabel('Timestep, k')
    ylabel('y(k)')
    grid on; hold off;
   subplot(2, 3, 2); hold on;
    acf(y, 25, 0.05, true, 25, true);
    title('ACF')
    xlabel('Lag')
    ylabel('ACF Coefficient')
    hold off;
subplot(2, 3, 3); hold on;
    pacf(y, 25, 0.05, true, true);
    title('PACF')
    xlabel('Lag')
    ylabel('PACF Coefficient')
    hold off;
subplot(2, 3, 4); hold on;
    scatter(y(1:N-1), y(2:N), '.')
    title('First Return Map')
    xlabel('y(k-1)')
    ylabel('y(k)')
    grid on; hold off;
subplot(2, 3, 5); hold on;
    ldf(y, 25, lambda, 4, 0.05, true);
    title('LDF')
    xlabel('Lag')
    ylabel('LDF Coefficient')
    hold off;
subplot(2, 3, 6); hold on;
    pldf(y, 25, lambda, 4, 0.05, true);
    title('PLDF')
    xlabel('Lag')
    ylabel('PLDF Coefficient')
    hold off;
sgtitle('Chaotic Logistic Map Test')

%% Test on White noise

N = 100;

y = randn(N, 1);

lambda = NaN; %0.1;

figure();
subplot(2, 3, 1); hold on;
    plot(y)
    title('Realization')
    xlabel('Timestep, k')
    ylabel('y(k)')
    grid on; hold off;
   subplot(2, 3, 2); hold on;
    acf(y, 25, 0.05, true, 25, true);
    title('ACF')
    xlabel('Lag')
    ylabel('ACF Coefficient')
    hold off;
subplot(2, 3, 3); hold on;
    pacf(y, 25, 0.05, true, true);
    title('PACF')
    xlabel('Lag')
    ylabel('PACF Coefficient')
    hold off;
subplot(2, 3, 4); hold on;
    scatter(y(1:N-1), y(2:N), '.')
    title('First Return Map')
    xlabel('y(k-1)')
    ylabel('y(k)')
    grid on; hold off;
subplot(2, 3, 5); hold on;
    ldf(y, 25, lambda, 4, 0.05, true);
    title('LDF')
    xlabel('Lag')
    ylabel('LDF Coefficient')
    hold off;
subplot(2, 3, 6); hold on;
    pldf(y, 25, lambda, 4, 0.05, true);
    title('PLDF')
    xlabel('Lag')
    ylabel('PLDF Coefficient')
    hold off;
sgtitle('White Noise Test')

plot_nl_residual(y, 'White Noise Test', 25, lambda, true, 0.05)