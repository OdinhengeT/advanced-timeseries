%% Non-Linear Timeseries Analysis Project 2023 - Implementation of Cox-de Boor's Recursive Algorithm with Derivatives
% Aron Paulsson and Torbj�rn Onshage
function [design_matrix] = cox_deBoor(x, knots, order, derivative_order)
    if nargin<3; order = 4; end
    if nargin<4; derivative_order = 0; end
    
    if ~isvector(x); error('cox_deBoor: input x must be vector'); end
    if ~isvector(knots); error('cox_deBoor: input knots must be vector'); end
    
    if ~iscolumn(x); x = x'; end
    if ~iscolumn(knots); knots = knots'; end
    
    if derivative_order > (order-1)
        design_matrix = zeros(length(x), length(knots) - order);
        return
    end
    
    Bs = cell(order, derivative_order+1);
    
    % First order B-splines 
    Bs{1,1} = zeros(length(x), length(knots) - 1);
    
    for i = 1:(length(knots) - 1)
        Bs{1,1}(:,i) = double( (x >= knots(i)) & (x < knots(i+1) | (knots(i+1) == knots(end))) );
        % The (knots(i+1) == knots(end)) is some necessary trickery for evaluation at the upper limit: max(knots)
    end

    % Cox-de Boor Recursive formula to calculate the requested order B-splines (with derivatives)
    for idx_p = 2:order
        for idx_k = 1:(derivative_order+1)
            if ~( (idx_k > idx_p) || ((derivative_order+1-idx_k) > (order - idx_p)) )
                Bs{idx_p,idx_k} = zeros(length(x), length(knots) - idx_p);
                
                for i = 1:(length(knots) - idx_p)
                    alpha = (x - knots(i)) ./ (knots(i+idx_p-1) - knots(i));
                    beta = (knots(i+idx_p) - x) ./ (knots(i+idx_p) - knots(i+1));
                    dAlpha = 1 ./ (knots(i+idx_p-1) - knots(i));
                    dBeta = -1 ./ (knots(i+idx_p) - knots(i+1));
                    
                    if isinf(dAlpha); alpha = 0.*x; dAlpha = 0; end
                    if isinf(dBeta); beta = 0.*x; dBeta = 0; end
                    
                    if idx_k == 1
                        Bs{idx_p,idx_k}(:,i) = alpha .* Bs{idx_p-1,idx_k}(:,i) + beta .* Bs{idx_p-1,idx_k}(:,i+1);
                        
                    elseif isempty( Bs{idx_p-1,idx_k} )
                        Bs{idx_p,idx_k}(:,i) = (idx_k-1) * dAlpha .* Bs{idx_p-1,idx_k-1}(:,i) ...
                                             + (idx_k-1) * dBeta .* Bs{idx_p-1,idx_k-1}(:,i+1);
                                         
                    else
                        Bs{idx_p,idx_k}(:,i) = (idx_k-1) * dAlpha .* Bs{idx_p-1,idx_k-1}(:,i) ...
                                             + alpha .* Bs{idx_p-1,idx_k}(:,i) ...
                                             + (idx_k-1) * dBeta .* Bs{idx_p-1,idx_k-1}(:,i+1) ...
                                             + beta .* Bs{idx_p-1,idx_k}(:,i+1); 
                    end
                end
            end
        end
    end
  
    design_matrix = Bs{end, end};
end