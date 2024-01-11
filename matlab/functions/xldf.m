%% Non-Linear Timeseries Analysis Project 2023 - Implementation of cross-LDF
% Aron Paulsson and Torbjörn Onshage
function [rho] = xldf(x, y, maxLag, lambda, order, signLvl, boolPlot, boolZeroLag)
    if nargin<3; maxLag = 10; end
    if nargin<4; lambda = NaN; end
    if nargin<5; order = 4; end
    if nargin<6; signLvl = 0.05; end
    if nargin<7; boolPlot = false; end
    if nargin<8; boolZeroLag = true; end
    
    if ~isvector(x); error('xldf: input x must be vector'); end
    if ~isvector(y); error('xldf: input y must be vector'); end
    
    if ~iscolumn(x); x = x'; end
    if ~iscolumn(y); y = y'; end
    
    if length(x) ~= length(y); error('xldf: inputs x and y must be same length'); end
    N = length(y);

    if (signLvl < 0) || (signLvl>1 )
        error('xldf: input signLvl must be between 0 and 1');
    else
        signScale = norminv( 1-signLvl/2, 0, 1 );
    end

    rho = ones(maxLag + (boolZeroLag), 1);

    SS0 = sum(( y(maxLag+1:N) - mean(y(maxLag+1:N)) ).^2);

    for lag = 1:maxLag

        [theta, knots] = smooth_spline( x(1+maxLag-lag:N-lag), y(1+maxLag:N), lambda, order );
        
        Bs = cell(length(knots), 1);
        for idx = 1:length(knots); Bs{idx} = cox_deBoor(x(1+maxLag-lag:N-lag), knots{idx}, order); end  
        
        residual = y(1+maxLag:N) - [Bs{:}] * theta;

        SS0k = sum( residual.^2 );
        
        R20 = (SS0-SS0k) / SS0;

        x_minmax = [min(x(1+maxLag-lag:N-lag)); max(x(1+maxLag-lag:N-lag))];
        for idx = 1:length(knots); Bs{idx} = cox_deBoor( x_minmax , knots{idx}, order); end  
        
        fk_minmax = [Bs{:}] * theta;
        
        rho(lag+1) = sign(fk_minmax(2)-fk_minmax(1))*sqrt((abs(R20)+R20)/2);
    end
    
    if boolPlot
        stem( int8(~boolZeroLag):maxLag, rho )
        xlabel('Lag')
        ylabel('Amplitude')
        condInt = signScale * ones(length(rho),1)/sqrt(length(y));
        condInt = condInt * sqrt( 1 + 2*sum( rho(1+int8(boolZeroLag):end).^2 ) );
        hold on
        plot( int8(~boolZeroLag):maxLag, condInt,'r--' )
        plot( int8(~boolZeroLag):maxLag, -condInt,'r--' )
        hold off
    end

end