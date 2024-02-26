function [fit2x1D,Y0Fit,Y1Fit,YAllFit]=MLE_Two1D_Gaussian_RealDataFixPSF(x,Y,est)
%%%This function fits two Gaussians to a line
% correctCoeffs = [ux1 ux2 sigmaValueX1 sigmaValueX2 b N1 N2];

if nargin<3
    
    close all
    SN=10;
    LL=100;
    x=1:LL;
    est.psf=LL/8;
    est.x0=LL/3;
    est.x1=2*LL/3;
    est.b0=0;
    est.N0=100;
    est.N1=50;
    Y0=est.N0*oneDGauss(x,est.x0,est.psf);
    Y1=est.N1*oneDGauss(x,est.x1,est.psf);
    Noise=1/SN*max(Y0+Y1)*randn(1,LL);
    Y=Y0+Y1+Noise;
end


x0 = [est.x0,est.x1,est.b0,est.N0,est.N1];

 % The function to be minimized is the negative of the log likelihood
datafun = @(params)((sum (expectedFixPSF (x,params,est.psf))) - sum (Y.*log (expectedFixPSF (x,params,est.psf))));
options = optimset ('MaxFunEvals', 10000, 'MaxIter', 10000, 'TolFun', 1e-5);
[paramsF,fval,exitflag,output]  = fminsearch (datafun, x0, options); %we make use of fminsearch to find the minimum of the -log-likelihood function
paramsF(4)=abs(paramsF(4));  %to ensure positive amplitudes (also in expectation expression)
paramsF(5)=abs(paramsF(5));  %to ensure positive amplitudes (also in expectation expression)

%Here we generate data with the fitted parameters.
Y0Fit =1/(est.psf*(2*pi)^0.5)*paramsF(4)*exp((-((x-paramsF(1)).^2./(2.*est.psf.^2)))) + paramsF(3); %Generate y-values with the fitted parameters
Y1Fit = 1/(est.psf*(2*pi)^0.5)*paramsF(5)*exp((-((x-paramsF(2)).^2./(2.*est.psf.^2))))+ paramsF(3); %Generate y-values with the fitted parameters
YAllFit = Y0Fit + Y1Fit - paramsF(3); %add the generated y-values to obtain the final curve

fit2x1D.x0=paramsF(1);
fit2x1D.x1=paramsF(2);
fit2x1D.B0=paramsF(3)
fit2x1D.N0=paramsF(4);  %between max 1 peak and two superimposed
fit2x1D.N1=paramsF(5);  %between max 1 peak and two superimposed
fit2x1D.psf=est.psf;     %these are not used

if nargin<3
    plot(x,Y0,'bo'); hold on;
    plot(x,Y0Fit,'b-'); hold on;
    plot(x,Y1,'ro');
    plot(x,Y1Fit,'r-');
    plot(x,Y,'ko');
    plot(x,YAllFit,'k-');
    legend('Peak1','Fit','Peak2','Fit','Combined', 'Sum Fits');
    
end



function p = oneDGauss (x,ux,s)
%This is the equation for a normalized1D gaussian
p = 1/(s*sqrt(2*pi))*exp (-(x-ux).^2./(2*s.^2));


function E = expectedFixPSF (x,params,psf)
Le=length(x);
E = abs(params (4)).*oneDGauss (x,params (1),psf)  +  abs(params (5)).*oneDGauss (x,params(2),psf) + abs(params (3));
