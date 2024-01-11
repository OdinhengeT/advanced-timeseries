clear all
close all
clc
%% Add subfolders
addpath('functions', 'data')
%%
dataTable_DK = readtable('Elspotprices_DK.csv');
data_strings_DK = dataTable_DK.SpotPriceDKK;
data_numeric = str2double(strrep(data_strings_DK, ',', '.'));
data_DK1 = data_numeric(1:2:end);
data_DK2 = data_numeric(2:2:end);
x_DK1 = 1:numel(data_DK1);
x_DK2 = 1:numel(data_DK2);

dataTable_SV = readtable('Elspotprices_SV.csv');
data_strings_SV = dataTable_SV.SpotPriceDKK;
data_numeric = str2double(strrep(data_strings_SV, ',', '.'));
data_SV3 = data_numeric(1:2:end);
data_SV4 = data_numeric(2:2:end);
x_SV3 = 1:numel(data_SV3);
x_SV4 = 1:numel(data_SV4);
%%
figure()
plot(x_DK1, data_DK1, 'r--');
hold on
plot(x_DK2, data_DK2, 'b--');
plot(x_SV3, data_SV3, 'g--');
plot(x_SV4, data_SV4, 'k--');
legend('DK1', 'DK2','SV3','SV4')
xlabel('Index');
ylabel('Data');
title('Plot of the Numeric Data');

%% 
N0 = 25000;
N = N0 + 3000;

y = data_SV4(N0:N);

figure();
subplot(2, 3, 1); hold on;
    plot(y)
    title('SE 4')
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
    scatter(y(1:end-1), y(2:end), '.')
    title('First Return Map')
    xlabel('y(k-1)')
    ylabel('y(k)')
    grid on; hold off;
subplot(2, 3, 5); hold on;
    ldf(y, 25, NaN, 4, 0.05, true);
    title('LDF')
    xlabel('Lag')
    ylabel('LDF Coefficient')
    hold off;
% subplot(2, 3, 6); hold on;
%     pldf(data_SV4, 25, NaN, 4, 0.05, true);
%     title('PLDF')
%     xlabel('Lag')
%     ylabel('PLDF Coefficient')
%     hold off;
    

%%

% model 1
Ap = [1 1 1 1 zeros(1, 19) 1 1 1];
Cp = [1 1 1 1 zeros(1, 20) 1];

model_init = idpoly(Ap, [], Cp);
model_init.Structure.a.Free = Ap;
model_init.Structure.c.Free = Cp;

model_pem = pem(iddata(y), model_init); 
present(model_pem); 

model_pem_res = resid(y, model_pem); model_pem_res = model_pem_res(length(Ap):end);
plot_residual(model_pem_res.y, 'PEM Model Residual', 50, false)

%% 
N = 5;
figure();
scatter(model_pem_res.y(1:end-N), model_pem_res.y(N+1:end), '.')
title([num2str(N), '-th Return Map'])
xlabel('y(k-2)')
ylabel('y(k)')
grid on; hold off;