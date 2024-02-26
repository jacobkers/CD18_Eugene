function [fit1x1D,YFit]=MLE_One1D_Gaussian_FreePSF(x,Y,est, plotit);

x0 = [est.x0,est.psf,est.b0,est.N0];
%The function to be minimized is the negative of the log likelihood
datafun = @(params)((sum (expected (x,params))) - sum (Y.*log (expected (x,params))));
options = optimset ('MaxFunEvals', 1000, 'MaxIter', 1000, 'TolFun', 1e-5);
[paramsF,fval,exitflag,output]  = fminsearch (datafun, x0, options); %we make use of fminsearch to find the minimum of the -log-likelihood function

%Here we generate data with the fitted parameters.
YFit = expected (x,paramsF); %Generate y-values with the fitted parameters

fit1x1D.x0=paramsF(1);
fit1x1D.N0=paramsF(4);
fit1x1D.b0=paramsF(3).^0.5;   %background content level
fit1x1D.psf=paramsF(2);



if plotit==1
    figure;
    plot(x,Y,'ok','LineWidth',1.5); hold on;
    plot(x,YFit,'r-','LineWidth',1.5);
    hold on;
    h = legend('Input Data','Fit of Data',5);
    set(h,'FontSize',14);
    xlabel ('Position ','fontsize',21);
    ylabel ('Intensity [a.u.]','fontsize',21);
    hold off;
    pause(0.2);
    close(gcf);
    [~]=ginput(1);
    end
end

function p = oneDGauss (x,ux,s)
%This is the equation for a 1D gaussian
p = 1/(sqrt(2*pi).*s) .* exp (-((x-ux).^2./(2*s.^2) ));
end


function E = expected (x,params)
% This is the expected value given the parameters
E = abs(params (4).*oneDGauss (x,params (1),params (2))) + params (3)^2 ;
end