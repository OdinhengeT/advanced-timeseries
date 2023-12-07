clear all
close all
clc
%% Add subfolders
addpath('data')
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