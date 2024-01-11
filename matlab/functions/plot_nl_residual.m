
function [] = plot_nl_residual(Residual, name, maxLag, lambda, zeroLag, signLvl)
if nargin < 2; name = '';                           end
if nargin < 3; maxLag = 25;                         end
if nargin < 4; lambda = NaN;                        end
if nargin < 5; zeroLag = true;                      end
if nargin < 6; signLvl = 0.05;                      end

if zeroLag; idx = 2; else; idx = 1; end

figure();
subplot(3, 2, [1 2])
    plot(Residual)
    title('Residual');
    grid on;
subplot(3, 2, 3)
    rho = ldf(Residual, maxLag, lambda, 4, signLvl, true, zeroLag);
    title('LDF')
subplot(3, 2, 4)
    phi = pldf(Residual, maxLag, lambda, 4, signLvl, true, zeroLag);
    title('PLDF')
subplot(3, 2, 5)
    normplot(rho(idx:end))
    if dagostinoK2test(rho(idx:end), signLvl);
        title('Normplot ACF, D''Agostino: True')
    else
        title('Normplot ACF, D''Agostino: False')
    end 
subplot(3, 2, 6)
    normplot(phi(idx:end))
    if dagostinoK2test(phi(idx:end), signLvl);
        title('Normplot PACF, D''Agostino: True')
    else
        title('Normplot PACF, D''Agostino: False')
    end 
sgtitle(name)




