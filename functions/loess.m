%% Non-Linear Timeseries Analysis Project 2023 - Implementation of LOESS
% Aron Paulsson and Torbjörn Onshage
function [y_eval, alpha] = loess(x, y, x_eval, alpha, order)
  
    if min(size(x)) ~= 1; error('LOESS: Input x must be vector'); end
    if min(size(y)) ~= 1; error('LOESS: Input y must be vector'); end 
    if ~iscolumn(x); x = x'; end
    if ~iscolumn(y); y = y'; end
    if length(x) ~= length(y)
        error('LOESS: Input x and y vectors must be of same length, are ' + str(length(x)) + ' and ' + str(length(y)))
    end
    
    if nargin<3; x_eval = NaN; end
    if nargin<4; alpha = NaN; end
    if nargin<5; order = 1; end
    
    N = length(x);
    
    if isnan(alpha)
        % Find "Optimal" alpha
        
        alphas = linspace((order+1)/N, 0.3, 5);
        MSEs = zeros(1, length(alphas));
        
        nbrTestPoints = round(max( 1, sqrt(0.1*N) + 5*(2/pi*atan(0.01*(N-150))+0.78) )); % Somewhat good number of points
        
        for j = 1:5
            idxTestPoints =  randperm(N, nbrTestPoints);

            % Exclude testpoints 
            x0 = x(setdiff(1:N, idxTestPoints));
            y0 = y(setdiff(1:N, idxTestPoints));

            for i = 1:length(alphas)
                [y_eval, ~] = loess( x0, y0, x(idxTestPoints), alphas(i) );
                MSEs(i) = MSEs(i) + sum( (y(idxTestPoints) - y_eval).^2 );
            end
        end
        clearvars x0 y0
        
        [~, idx] = min(MSEs);
        alpha = alphas(idx);
    end
    
    if ~isnan(x_eval)
        % Evaluate at specified point
        
        nbrFitPoints = 1 + order + round(alpha * N);
        
        y_eval = zeros(length(x_eval), 1);
        
        for i = 1:length(x_eval) 

            [~, idxs] = minVals((x - x_eval(i)).^2, nbrFitPoints+1);
            
            dist_max = abs(x( idxs(end) ) - x_eval(i));
            idxs = idxs(1:(end-1));
            
            X = ones(nbrFitPoints, order+1);
            for j = 1:order
               X(:, j+1) = x(idxs).^j;
            end
            
            W = diag( kernel( (x(idxs)-x_eval(i)) / dist_max ) );
 
            beta = (X'*W*X) \ X'*W*y(idxs);
            
            X_eval = ones(1, order+1);
            for j = 1:order
               X_eval(:, j+1) = x_eval(i).^j;
            end
            
            y_eval(i) = X_eval*beta;
        end
    end
end

function [minVals, minIdxs] = minVals(y, N)
     [y_sort, y_idx] = sort(y);
     minVals = y_sort(1:N);
     minIdxs = y_idx(1:N);
end

function K = kernel(u)
    K = 0.75.*(1 - u.^2);
end
