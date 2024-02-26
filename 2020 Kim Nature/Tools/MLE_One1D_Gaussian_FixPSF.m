function [fit1x1Dfix,YFit]=MLE_One1D_Gaussian_FixPSF(x,Y,est,plotit);
%JWJK_B:-------------------------------------------------------------------
%Title: MLE_One1D_Gaussian_FixPSF; 
%
%Summary: %This function fits one 1D Gaussian to a line, fixed width
%correctCoeffs = [ux1 ux2 sigmaValueX1 sigmaValueX2 b N1 N2];
%
%:JWJK_B-------------------------------------------------------------------  


a=1; % this should in principle be the pixel size
%Defining the initial conditions for the optimization procedure
x0 = [est.x0,est.b0,est.N0];


 % The function to be minimized is the negative of the log likelihood
datafun = @(params)((sum (expectedFixPSF (x,params,est.psf))) - sum (Y.*log (expectedFixPSF (x,params,est.psf))));
options = optimset ('MaxFunEvals', 10000, 'MaxIter', 10000, 'TolFun', 1e-5);
[paramsF,fval,exitflag,output]  = fminsearch (datafun, x0, options); %we make use of fminsearch to find the minimum of the -log-likelihood function
paramsF(3)=abs(paramsF(3));  %to ensure positive amplitudes (also in expectation expression)

%Here we generate data with the fitted parameters.
YFit = expectedFixPSF (x,paramsF,est.psf); %Generate y-values with the fitted parameters


fit1x1Dfix.x0=paramsF(1);
fit1x1Dfix.b0=paramsF(2).^0.5;   %background content
fit1x1Dfix.N0=paramsF(3);
fit1x1Dfix.psf=est.psf;

if plotit==1
    figure;
    plot(x,Y,'o-k','LineWidth',1.5); hold on;
    plot(x,YFit,'o-r','LineWidth',1.5);
    hold on;
    h = legend('Input Data','Fit of Data');
    set(h,'FontSize',14);
    xlabel ('Position ','fontsize',21);
    ylabel ('Intensity [a.u.]','fontsize',21);
    hold off;
    pause(0.2);
    close(gcf);
    end
end

function p = oneDGauss (x,ux,s)
%This is the equation for a normalized1D gaussian
p = 1/(s*sqrt(2*pi))*exp (-(x-ux).^2./(2*s.^2));
end


function E = expectedFixPSF (x,params,psf)
% This is the expected value given the parameters
E = abs(params (3)).*oneDGauss (x,params (1),psf) + abs(params (2)).^2;
end