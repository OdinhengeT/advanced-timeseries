%% Non-Linear Timeseries Analysis Project 2023 - Implementation of Smoothing Splines
% Aron Paulsson and Torbj�rn Onshage
function [theta, knots, lambda] = smooth_spline(x, y, lambda, order, x_range)
    if nargin<3; lambda = NaN; end
    if nargin<4; order = 4; end
    if nargin<5; x_range = [min(x); max(x)] + 0.1.* [-1; 1] * (max(x) - min(x)); end
    
    if ~isvector(y); error('smooth_spline: input y must be vector'); end
    N = length(y); % Number of samples
    
    if max(size(x)) ~= N; error('smooth_spline: input x must be of same length as input y'); end
    
    if size(x, 1) ~= N; x = x'; end
    if ~iscolumn(y); y = y'; end
    
    if size(x_range, 1) ~= 2 || size(x_range, 2) ~= size(x, 2)
        error('smooth_spline: input x_range must have 2 rows and the same number of columns as x'); 
    end
    
    if ~isvector(order); error("smooth_spline: input order must be vector or scalar"); end
    
    if isnan(order); order = 4.* ones(size(x, 2), 1); end
        
    if length(order) ~= size(x, 2)
        if length(order) == 1
            order = order .* ones(size(x, 2), 1); 
        else  
            error("smooth_spline: input order must be scalar or same length as the number of inputs (columns) in x")
        end
    end
    
    if isnan(lambda)
        
        lambdas = [0.1, 0.01, 0.001, 0.0001, 0.00001];
        MSEs = zeros(1, length(lambdas));
        
        nbrTestPoints = ceil(min([ sqrt(0.3.*N), 200 + N.^0.2 ]));

        for repeats = 1:5
            idxTestPoints =  randperm(N, nbrTestPoints);

            % Exclude testpoints 
            x_model = x(setdiff(1:N, idxTestPoints), :);
            y_model = y(setdiff(1:N, idxTestPoints), :);
            
            x_test = x(union(1:N, idxTestPoints), :);
            y_test = y(union(1:N, idxTestPoints), :);
            
            for idx_lambda = 1:length(lambdas)
            
                [theta, knots] = smooth_spline(x_model, y_model, lambdas(idx_lambda), order, x_range);
                
                Bs = cell(length(knots), 1);
                for idx = 1:length(knots); Bs{idx} = cox_deBoor(x_test(:, idx), knots{idx}, order(idx)); end  
                
                MSEs(idx_lambda) = MSEs(idx_lambda) + sum( (y_test - [Bs{:}] * theta).^2 );
            end
        end
        
        [~, idx] = min(MSEs);
        lambda = lambdas(idx);
    end
    
    knots = cell(size(x, 2), 1);
    omegas = cell(size(x, 2), 1);
    Bs = cell(size(x, 2), 1);
    
    for idx_x = 1:size(x, 2)
        K = ceil( min([ N/3, 3*sqrt(N), 200 + 2.*N.^0.2 ]) ); % Number of internal knots

        internal_knots = unique(quantile(x(:,idx_x), (1:K)./(K+1))');
        knots{idx_x} = [x_range(1, idx_x).*ones(order(idx_x),1); internal_knots; x_range(2, idx_x).*ones(order(idx_x),1)];

        omegas{idx_x} = penalty_matrix(knots{idx_x}, order(idx_x));

        Bs{idx_x} = cox_deBoor(x(:, idx_x), knots{idx_x}, order(idx_x));
    end
    
    omega = blkdiag(omegas{:});
    B = [Bs{:}];
    
    theta = (B'*B + N.*lambda.*omega) \ B' * y;
end

% Calculates the Penalty matrix introduced by O'Sullivan,
% but uses Gauss-Legendre Quadrature instead of Simpson's rule
function [omega] = penalty_matrix(knots, order)
    if nargin<2; order = 4; end
    
    if ~isvector(knots); error('penalty_matrix: input knots must be vector'); end
    if ~iscolumn(knots); knots = knots'; end
    
    % Calculate nodes and weights for Gauss-Legendre quadrature numerical integration
    [nodes, weights] = lgwt( max(1, order-2), 0, 1);
    
    nodes = nodes(end:-1:1);
    weights = weights(end:-1:1);
    
    dKnots = knots(2:end) - knots(1:end-1);
    
    x_tilde = repelem(knots(1:end-1), length(nodes), 1) + reshape( nodes .* dKnots', [] , 1);
    
    weights = reshape( weights .* dKnots', [] , 1);
    
    B = cox_deBoor(x_tilde, knots, order, 2);
    
    omega = B'*diag(weights)*B;
end

% Enormous Thanks to Greg von Winckel for publishing this script
% https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes
% 
% Copyright (c) 2009, Greg von Winckel
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

function [x,w] = lgwt(N,a,b)
    % lgwt.m
    %
    % This script is for computing definite integrals using Legendre-Gauss 
    % Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
    % [a,b] with truncation order N
    %
    % Suppose you have a continuous function f(x) which is defined on [a,b]
    % which you can evaluate at any x in [a,b]. Simply evaluate it at all of
    % the values contained in the x vector to obtain a vector f. Then compute
    % the definite integral using sum(f.*w);
    %
    % Written by Greg von Winckel - 02/25/2004
    N=N-1;
    N1=N+1; N2=N+2;
    xu=linspace(-1,1,N1)';
    % Initial guess
    y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
    % Legendre-Gauss Vandermonde Matrix
    L=zeros(N1,N2);
    % Derivative of LGVM
    Lp=zeros(N1,N2);
    % Compute the zeros of the N+1 Legendre Polynomial
    % using the recursion relation and the Newton-Raphson method
    y0=2;
    % Iterate until new points are uniformly within epsilon of old points
    while max(abs(y-y0))>eps
        L(:,1)=1;
        Lp(:,1)=0;

        L(:,2)=y;
        Lp(:,2)=1;

        for k=2:N1
            L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
        end

        Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   

        y0=y;
        y=y0-L(:,N2)./Lp;

    end
    % Linear map from[-1,1] to [a,b]
    x=(a*(1-y)+b*(1+y))/2;      
    % Compute the weights
    w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;

end