%% Non-Linear Timeseries Analysis Project 2023 - Implementation of cross-PLDF
% Aron Paulsson and Torbjörn Onshage
function [phi] = xpldf(x, y, maxLag, lambda, order, signLvl, boolPlot, boolZeroLag)
    if nargin<3; maxLag = 10; end
    if nargin<4; lambda = NaN; end
    if nargin<5; order = 4; end
    if nargin<6; signLvl = 0.05; end
    if nargin<7; boolPlot = false; end
    if nargin<8; boolZeroLag = true; end
    
    if ~isvector(x); error('xpldf: input x must be vector'); end
    if ~isvector(y); error('xpldf: input y must be vector'); end
    
    if ~iscolumn(x); x = x'; end
    if ~iscolumn(y); y = y'; end
    
    if length(x) ~= length(y); error('xpldf: inputs x and y must be same length'); end
    N = length(y);

    if (signLvl < 0) || (signLvl>1 )
        error('xpldf: input signLvl must be between 0 and 1');
    else
        signScale = norminv( 1-signLvl/2, 0, 1 );
    end

    phi = ones(maxLag + (boolZeroLag), 1);

    SS0ks = zeros(maxLag+1);    
    SS0ks(1) = sum(( y(maxLag+1:N) - mean(y(maxLag+1:N)) ).^2); % SS0
    
    X = [];
    
    for lag = 1:maxLag
        
        X = [X, x(1+maxLag-lag:N-lag)];
        [theta, knots] = smooth_spline(X, y(1+maxLag:N), lambda, order );
        
        Bs = cell(length(knots), 1);
        for idx = 1:length(knots); Bs{idx} = cox_deBoor(X(:, idx), knots{idx}, order); end  
        
        residual = y(1+maxLag:N) - [Bs{:}] * theta;
        
        SS0ks(1+lag) = sum( residual.^2 );
        
        R20 = (SS0ks(lag) - SS0ks(1+lag)) / SS0ks(lag);

        fkk_minmax = [min(x(1+maxLag-lag:N-lag)); max(x(1+maxLag-lag:N-lag))];
        
        B = cox_deBoor( fkk_minmax , knots{end}, order );
        
        yhat_minmax = B * theta( end-size(B,2)+1:end );
        
        phi(lag+1) = sign(yhat_minmax(2)-yhat_minmax(1))*sqrt((abs(R20)+R20)/2);
    end
    
    if boolPlot
        stem( int8(~boolZeroLag):maxLag, phi )
        xlabel('Lag')
        ylabel('Amplitude')
        condInt = signScale * ones(length(phi),1)/sqrt(length(y));
        condInt = condInt * sqrt( 1 + 2*sum( phi(1+int8(boolZeroLag):end).^2 ) );
        hold on
        plot( int8(~boolZeroLag):maxLag, condInt,'r--' )
        plot( int8(~boolZeroLag):maxLag, -condInt,'r--' )
        hold off
    end

end

