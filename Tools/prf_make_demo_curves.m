 function prf=prf_make_demo_curves(whatcurve);
 %make fake curves for demo options
if nargin<1,whatcurve= 'multipeaks';  end
switch whatcurve
    case 'plectoneme'
        axz=1:100; prf=0*axz;
        for ii=20:10:80
            prf=prf+exp(-((axz-ii)/10).^2);
            prf=prf+0.2*exp(-((axz-50)/5).^2);
        end  
    case 'multipeaks'
        axz=1:200; sig=1; prf=randn(1,200);
        mus=rand(6,1)*160+20;
        for ii=1:6
            mu=mus(ii); 
            prf=prf+10*exp(-((axz-mu)/(2*sig)).^2);
        end  
end
if nargin<1,close all; plot(prf);  end
