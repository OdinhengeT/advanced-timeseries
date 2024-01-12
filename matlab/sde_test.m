%% SDE 1
clear all;
close all;
addpath('functions', 'data');
clc;

A = 100;
c = 5;
w = 6.29;
mu = 0.1;

param = [A, c, w, mu];

func_derivative = @(t, X, param) [
    -param(3) .* X(2,:) - param(4) .* X(1,:) .* ( X(1,:).^2 + X(2,:).^2 - 1 );
     param(3) .* X(1,:) - param(4) .* X(2,:) .* ( X(1,:).^2 + X(2,:).^2 - 1 )
];

dt = 1/12000;

transitionFunction = @(t, X, param, dt) EEstep(func_derivative, param, 0, X, dt);

observationFunction = @(t, X, param) param(1)*exp(-2.*param(2)) .* ( exp( param(2) .* (1+X(2,:)) ) - 1 );

T = 25;

X = zeros( 2, ceil(T/dt)+1 );
X(:,1) = [-1, 0];

Y = zeros( 1, ceil(T/dt)+1 );
Y(:,1) = observationFunction(0, X(:,1), param);

for t = 2:length(Y)
    
    X(:,t) = transitionFunction(0, X(:,t-1), param, dt) + 0.025 .* sqrt(dt) .* randn(2,1);
    Y(:,t) = observationFunction(0, X(:,t-1), param) + 0.25 .* randn(1);
   
end

t = 0:dt:T;

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on;
plot(t, X(1,:))
plot(t, X(2,:))
title('Realizations')
legend('X1', 'X2')
xlabel('time')
ylabel('Value')
grid on; hold off;
subplot(1,2,2); hold on;
plot(X(1,:), X(2,:))
title('Phase Portrait')
xlabel('X1')
ylabel('X2')
grid on; hold off;
sgtitle({'Simulation of SDE 1'})

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on;
plot(t, Y)
title('Realizations')
legend('Y')
xlabel('time')
ylabel('Value')
grid on; hold off;
subplot(1,2,2); hold on;
plot(X(1,:), Y)
title('Phase Portrait')
xlabel('X1')
ylabel('Y')
grid on; hold off;
sgtitle({'Simulation of SDE 1'})

%% SDE 2
clear all;
close all;
addpath('functions', 'data');
clc;

A = 100;
c = 5;
w = 6.29;
mu = 0.1;

param = [A, c, w, mu];

func_derivative = @(t, X, param) [
    -param(3) .* X(2,:) - param(4) .* X(1,:) .* ( X(1,:).^2 + X(2,:).^2 - 1 );
     param(3) .* X(1,:) - param(4) .* X(2,:) .* ( X(1,:).^2 + X(2,:).^2 - 1 );
     param(2) .* param(3) .* X(1,:) .* exp( param(2) .* (1+X(2,:)) ) .* A .* exp(-2.*param(2))
];

dt = 1/120;

transitionFunction = @(t, X, param, dt) RK4step(func_derivative, param, 0, X, dt);

observationFunction = @(t, X, param) X(3,:);

T = 25;

X = zeros( 3, ceil(T/dt)+1 );
X(:,1) = [-1, 0, 0.02];

Y = zeros( 1, ceil(T/dt)+1 );
Y(:,1) = observationFunction(0, X(:,1), param);

for t = 2:length(Y)
    
    X(:,t) = transitionFunction(0, X(:,t-1), param, dt) + 0.025 .* sqrt(dt) .* randn(3,1);
    Y(:,t) = observationFunction(0, X(:,t-1), param) + 0.25 .* randn(1);
   
end

t = 0:dt:T;

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on;
plot(t, X(1,:))
plot(t, X(2,:))
title('Realizations')
legend('X1', 'X2')
xlabel('time')
ylabel('Value')
grid on; hold off;
subplot(1,2,2); hold on;
plot(X(1,:), X(2,:))
title('Phase Portrait')
xlabel('X1')
ylabel('X2')
grid on; hold off;
sgtitle({'Simulation of SDE 2'})

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on;
plot(t, Y)
title('Realizations')
legend('Y')
xlabel('time')
ylabel('Value')
grid on; hold off;
subplot(1,2,2); hold on;
plot(X(1,:), Y)
title('Phase Portrait')
xlabel('X1')
ylabel('Y')
grid on; hold off;
sgtitle({'Simulation of SDE 2'})


%% Other Testing
clear all;
close all;
clc;

param = [3.14, 0.45, 6.29];

func_derivative = @(t, X, param) [
    param(2) .* ( param(1) - X(1,:) );
    param(3);
];

dt = 1/12;

T = 30;
Y_0 = [3.5, -1.8, -1];

Y = zeros( 3, ceil(T/dt)+1 );
Y(:,1) = Y_0';

sigma_r = exp(-6.4);
sigma_theta = exp(-1.4);
sigma_Z = exp(1.75);

sigma = [sigma_r, sigma_theta, sigma_Z]';

for t = 2:length(Y)
    Y(:,t) = RK4step(func_derivative, param, t*dt, Y(:,t-1), dt) + sqrt(dt) .* sigma .* randn(3,1);
end

t = 0:dt:T;

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on;
plot(t, Y(1,:) )
plot(t, mod(Y(2,:), 2*pi) )
plot(t, Y(3,:))
title('Realizations')
legend('R', 'w', 'Z')
xlabel('time')
ylabel('Value')
grid on; hold off;
subplot(1,2,2); hold on;
plot( Y(1,:).*cos(Y(2,:)), Y(1,:).*sin(Y(2,:)) )
title('Phase Portrait')
xlabel('Y1')
ylabel('Y2')
grid on; hold off;
sgtitle({'Simulation of my SDE'})

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on;
plot(t, 0.75 .* (exp(Y(1,:) .* (1 + sin(Y(2,:)))) - 1))
title('Realizations')
legend('Y3')
xlabel('time')
ylabel('Value')
grid on; hold off;
subplot(1,2,2); hold on;
plot(Y(1,:).*cos(Y(2,:)), 0.75 .* (exp(Y(1,:) .* (1 + sin(Y(2,:)))) - 1))
title('Phase Portrait')
xlabel('Y1')
ylabel('Y2')
grid on; hold off;
sgtitle({'Simulation of my SDE'})


%% Another attempt
clear all;
close all;
clc;

A = 3;
c = 5;
w = 6.25;

param = [A, c, w];

func_derivative = @(t, X, param) [
    -param(3).*X(:,2) , param(3).*X(:,1) , param(2).*param(3).*X(:,1) .* exp( param(2).*( 1 + X(:,2) ) ) .* param(1).*exp(-2.*param(2));
];

dt = 1/12;

T = 10;
X_0 = [0.0, -1.0, 0.02];

X = zeros( ceil(T/dt)+1 , 3 );
X(1,:) = X_0;

sigma = 0.01;

for t = 2:length(X)
    %X(:,t) = X(:,t-1) + dt * func_derivative((t-1)*dt, Y(:,t-1), param) + sqrt(dt) .* sigma .* randn(1,1);
    
    %param(4) = param(4) +  sqrt(dt) .* sigma .* randn(1,1);
    
    X(t, :) = RK4step(func_derivative, param, (t-1)*dt, X(t-1, :), dt) + sqrt(dt) .* sigma .* randn(1,3);
    
%      if X(t, 3) < 0
%          X(t,3) = 0;
%      end
end

t = 0:dt:T;

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on;
plot(t, X(:,1) )
plot(t, X(:,2) )
plot(t, X(:,3))
title('Realizations')
legend('X1', 'X2', 'X3')
xlabel('time')
ylabel('Value')
grid on; hold off;
subplot(1,2,2); hold on;
plot( X(:,1), X(:,2) )
title('Phase Portrait')
xlabel('X1')
ylabel('X2')
grid on; hold off;
sgtitle({'Simulation of my SDE'})

figure('Position', [100, 100, 800, 400]);
subplot(1,2,1); hold on;
plot(t, X(:,3))
title('Realizations')
legend('X3')
xlabel('time')
ylabel('Value')
grid on; hold off;
subplot(1,2,2); hold on;
plot(X(:,1), X(:,3))
title('Phase Portrait')
xlabel('X1')
ylabel('X3')
grid on; hold off;
sgtitle({'Simulation of my SDE'})
