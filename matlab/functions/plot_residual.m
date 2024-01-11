
function [] = plot_residual(Residual, name, maxLag, zeroLag, signLvl)
if nargin < 5; signLvl = 0.05;                      end
if nargin < 4; zeroLag = true;                      end
if nargin < 3; maxLag = round(length(Residual)/4);  end
if nargin < 2; name = '';                           end

if zeroLag; idx = 2; else; idx = 1; end

figure();
subplot(3, 2, [1 2])
    plot(Residual)
    title('Residual');
    grid on;
subplot(3, 2, 3)
    rho = acf(Residual, maxLag, signLvl, 1, 0, zeroLag);
    title('ACF')
subplot(3, 2, 4)
    phi = pacf(Residual, maxLag, signLvl, 1, zeroLag);
    title('PACF')
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




