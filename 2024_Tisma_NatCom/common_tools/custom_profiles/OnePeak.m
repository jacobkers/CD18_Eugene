function PP= OnePeak(x,ux,s,PeakMode)
%JWJK_C:-------------------------------------------------------------------
%Title: One Gaussian Peak
%
%Summary: %This is the equation for a normalized1D gaussian peak value one,
%optional periodic boundaries

%Project:JacobKers 
%:JWJK_C-------------------------------------------------------------------
LL=length(x);
switch PeakMode
    case 'Nonperiodic'
        PP =exp (-(x-ux).^2./(2*s.^2));
    case 'Periodic'
        PP= exp (-(x-(ux-LL)).^2./(2*s.^2))+...
            exp (-(x-ux).^2./(2*s.^2))+...
            exp (-(x-(ux+LL)).^2./(2*s.^2));
end