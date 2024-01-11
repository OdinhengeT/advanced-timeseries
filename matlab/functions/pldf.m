%% Non-Linear Timeseries Analysis Project 2023 - Implementation of PLDF
% Aron Paulsson and Torbjörn Onshage
function [phi] = pldf(y, maxLag, lambda, order, signLvl, boolPlot, boolZeroLag)
    if nargin<2; maxLag = 10; end
    if nargin<3; lambda = NaN; end
    if nargin<4; order = 4; end
    if nargin<5; signLvl = 0.05; end
    if nargin<6; boolPlot = false; end
    if nargin<7; boolZeroLag = true; end
    
    phi = xpldf(y, y, maxLag, lambda, order, signLvl, boolPlot, boolZeroLag);
end

