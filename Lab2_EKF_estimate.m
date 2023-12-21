%% Part 4
clear all
close all
clc

%% Simulate data
a = 0.4;
N = 10000;
num_simulations = 20;
Y = zeros(num_simulations, N);
extra_n = 100;

for i = 1:num_simulations
    v = randn(1, N + extra_n);
    e = randn(1, N + extra_n);
    x = zeros(1, N + extra_n);
    y = zeros(1, N + extra_n);

    x(1) = v(1);
    y(1) = x(1) + e(1);

    for t = 2:N + extra_n
        x(t) = a * x(t-1) + v(t);
        y(t) = x(t) + e(t);
    end
    y = y(extra_n+1:end);
    Y(i,:) = y; % Store each simulation in a row of the matrix X
end
Y = Y';

plot(Y)
%%
numParticles = [1000 2];
ObservationDensity = @(y,state) normpdf(y,state(:,1),1); 
TransitionDensity = @(state, t) [state(:,2).*state(:,1), state(:,2)] + [1, 1/(t^(2/3))] .* randn(numParticles);  % Transition density,   q(x_k+1 | x_k)

AEst = zeros(N,20);
for i = 1:20 
    state = SMC_estimate(Y(:,i), ObservationDensity, TransitionDensity, numParticles);
    AEst(:,i) = state(:,2);
end
figure()
plot(AEst)
xlabel("Iteration",'fontweight','bold','fontsize',12);
ylabel("â",'fontweight','bold','fontsize',12);
title("Estimation of the parameter a = 0.4");
grid on



%%
% aInit = 0.5;
% aVarInit = 1;
% sigma_e = 1;
% sigma_v = 10;
% 
% f = @(x)[(x(2)*x(1)), x(2)];
% h = @(x)(x(1));
% figure()
% hold on
% for k = 1:20
%     obj = unscentedKalmanFilter(f,h,[0,aInit]);
%     obj.ProcessNoise = [sigma_v, 0; 0, 0];
%     obj.MeasurementNoise = sigma_e;
%     obj.StateCovariance = [sigma_e, 0; 0, 10];
% 
%     x = zeros(2,N);
%     for i = 1:N
%         x(:,i) = correct(obj,Y(i,k));
%         predict(obj);
%     end
%     ahat = x(2,:);
%     plot(ahat)
% end
% xlabel("Iteration",'fontweight','bold','fontsize',12)
% ylabel("â",'fontweight','bold','fontsize',12)
% title(["Estimation of the parameter a = 0.4 with aInit = 0.5,", "\sigma_v^2 = 10 and inital a variance = 1"])    %, for 20 different", "realizations of a model with \sigma_v^2 = 10 initial V(â) = 1"])
% grid on

%% a = 0.5, sigma_v = 10, aVarInit = 1
aInit = 0.5;
aVarInit = 1;
sigma_v = 10;

AEst = zeros(N,20);
AVarEst = zeros(N,20);
for i = 1:20 
    [y, avar] = ekf_estimate(Y(:,i), aInit, aVarInit, sigma_v);
    AEst(:,i) = y(:,2);
    AVarEst(:,i) = avar;
end
figure()
plot(AEst)
xlabel("Iteration",'fontweight','bold','fontsize',12)
ylabel("â",'fontweight','bold','fontsize',12)
title(["Estimation of the parameter a = 0.4 with aInit = 0.5,", "\sigma_v^2 = 10 and inital a variance = 1"])    %, for 20 different", "realizations of a model with \sigma_v^2 = 10 initial V(â) = 1"])
grid on

%figure()
%plot(AVarEst);

%% a = -0.5, sigma_v = 10, aVarInit = 1
aInit = -0.5;
aVarInit = 1;
sigma_v = 10;

AEst = zeros(N,20);
for i = 1:20 
    [y, avar] = ekf_estimate(Y(:,i), aInit, aVarInit, sigma_v);
    AEst(:,i) = y(:,2);
end
figure()
plot(AEst)
xlabel("Iteration",'fontweight','bold','fontsize',12)
ylabel("â",'fontweight','bold','fontsize',12)
title(["Estimation of the parameter a = 0.4 with aInit = -0.5,", " \sigma_v^2 = 10 and inital a variance = 1"])
grid on

%% a = 0.5, sigma_v = 1, aVarInit = 1
aInit = 0.5;
aVarInit = 1;
sigma_v = 1;

AEst = zeros(N,20);
for i = 1:20 
    [y, avar] = ekf_estimate(Y(:,i), aInit, aVarInit, sigma_v);
    AEst(:,i) = y(:,2);
end
figure()
plot(AEst)
xlabel("Iteration",'fontweight','bold','fontsize',12)
ylabel("â",'fontweight','bold','fontsize',12)
title(["Estimation of the parameter a = 0.4 with aInit = 0.5,", " \sigma_v^2 = 1 and inital a variance = 1"])
grid on

%% a = -0.5, sigma_v = 1, aVarInit = 1
aInit = -0.5;
aVarInit = 1;
sigma_v = 1;

AEst = zeros(N,20);
for i = 1:20 
    [y, avar] = ekf_estimate(Y(:,i), aInit, aVarInit, sigma_v);
    AEst(:,i) = y(:,2);
end
figure()
plot(AEst)
xlabel("Iteration",'fontweight','bold','fontsize',12)
ylabel("â",'fontweight','bold','fontsize',12)
title(["Estimation of the parameter a = 0.4 with aInit = -0.5,", " \sigma_v^2 = 1 and inital a variance = 1"])
grid on

%% a = 0.5, sigma_v = 10, aVarInit = 10
aInit = 0.5;
aVarInit = 10;
sigma_v = 10;

AEst = zeros(N,20);
for i = 1:20 
    [y, avar] = ekf_estimate(Y(:,i), aInit, aVarInit, sigma_v);
    AEst(:,i) = y(:,2);
end
figure()
plot(AEst)
xlabel("Iteration",'fontweight','bold','fontsize',12)
ylabel("â",'fontweight','bold','fontsize',12)
title(["Estimation of the parameter a = 0.4 with aInit = 0.5,", " \sigma_v^2 = 10 and inital a variance = 10"])
grid on

%% a = -0.5, sigma_v = 10, aVarInit = 10
aInit = -0.5;
aVarInit = 10;
sigma_v = 10;

AEst = zeros(N,20);
for i = 1:20 
    [y, avar] = ekf_estimate(Y(:,i), aInit, aVarInit, sigma_v);
    AEst(:,i) = y(:,2);
end
figure()
plot(AEst)
xlabel("Iteration",'fontweight','bold','fontsize',12)
ylabel("â",'fontweight','bold','fontsize',12)
title(["Estimation of the parameter a = 0.4 with aInit = -0.5,", " \sigma_v^2 = 10 and inital a variance = 10"])
grid on

%% a = 0.5, sigma_v = 1, aVarInit = 10
aInit = 0.5;
aVarInit = 10;
sigma_v = 1;

AEst = zeros(N,20);
for i = 1:20 
    [y, avar] = ekf_estimate(Y(:,i), aInit, aVarInit, sigma_v);
    AEst(:,i) = y(:,2);
end
figure()
plot(AEst)
xlabel("Iteration",'fontweight','bold','fontsize',12)
ylabel("â",'fontweight','bold','fontsize',12)
title(["Estimation of the parameter a = 0.4 with aInit = 0.5,", " \sigma_v^2 = 1 and inital a variance = 10"])
grid on

%% a = -0.5, sigma_v = 1, aVarInit = 10
a = -0.5;
aVarInit = 10;
sigma_v = 1;

AEst = zeros(N,20);
for i = 1:20 
    [y, avar] = ekf_estimate(Y(:,i), a, aVarInit, sigma_v);
    AEst(:,i) = y(:,2);
end
figure()
plot(AEst)
xlabel("Iteration",'fontweight','bold','fontsize',12)
ylabel("â",'fontweight','bold','fontsize',12)
title(["Estimation of the parameter a = 0.4 with aInit = -0.5,", " \sigma_v^2 = 1 and inital a variance = 10"])
grid on