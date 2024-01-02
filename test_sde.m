clear all
close all
addpath('functions', 'data');
clc

delta_t = 0.0001;
T = 10;

sigma = 0.00001;

x0 = 3;
x = zeros( 1, ceil(T/delta_t)+1 );

x(:, 1) = x0;

a = -0.3;

for t = 2:length(x)
    x(:, t) = (1 + a * delta_t) * x(:,t-1) + sigma * sqrt(delta_t) * randn(1);
end

t = 0:delta_t:T;

figure(); hold on;
plot(t, x(1,:), '-')
plot(t, x0.*exp(a.*t), 'k--')
hold off;

figure(); hold on;
scatter(x(1, 1:end-1), x(1, 2:end))

dx =  (x(1, 3:end) -  x(1, 1:end-2)) / (2 .* delta_t);

figure(); hold on;
scatter(x(1, 2:end-1), dx)

%%
clear all;
close all;
clc;

param = [3.2, 0.02, 6.29];

func_derivative = @(t, X, param) [
    -param(3) .* X(2,:) - 0.05 .* X(1,:) .* ( X(1,:).^2 + X(2,:).^2 - param(1).^2 );
     param(3) .* X(1,:) - 0.05 .* X(2,:) .* ( X(1,:).^2 + X(2,:).^2 - param(1).^2 )
];

dt = 1/12;

T = 30;
Y_0 = [param(1), 0];

Y = zeros( 2, ceil(T/dt)+1 );
Y(:,1) = Y_0;

sigma = 0.2;

for t = 2:length(Y)
    
    rad_noise = 0.1 .* sqrt(dt) .* randn(1) .* Y(:, t-1)./norm(Y(:, t-1));
    
    phase_noise = 0.2 .* sqrt(dt) .* randn(1) .* [Y(2, t-1); -Y(1, t-1)]./norm(Y(:, t-1));
    
    Y(:,t) = RK4step(func_derivative, param, t*dt, Y(:,t-1), dt) + rad_noise + phase_noise;
end

t = 0:dt:T;
Y3 = param(2).* ( exp(param(1) + Y(2,:) ) - 1 );

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on;
plot(t, Y(1,:))
plot(t, Y(2,:))
plot(t, param(1)*cos(param(3)*t), '--')
title('Realizations')
legend('Y1', 'Y2', 'a cos(wt)')
xlabel('time')
ylabel('Value')
grid on; hold off;
subplot(1,2,2); hold on;
plot(Y(1,:), Y(2,:))
title('Phase Portrait')
xlabel('Y1')
ylabel('Y2')
grid on; hold off;
sgtitle({'Simulation of my SDE'})

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on;
plot(t, Y3)
title('Realizations')
legend('Y3')
xlabel('time')
ylabel('Value')
grid on; hold off;
subplot(1,2,2); hold on;
plot(Y(1,:), Y3)
title('Phase Portrait')
xlabel('Y1')
ylabel('Y2')
grid on; hold off;
sgtitle({'Simulation of my SDE'})