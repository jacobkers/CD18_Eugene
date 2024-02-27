function msd=Get_MeanSquareDisplacement(x,y,lag);
%JWJK_C:
%Mean Square Displacment.
%
% %This function calculates the mean square displacement of a xy-trajectory
% %afo lag time in frames.
%:JWJK_C
lx=length(x);
if lag<lx-1                     %only calculate msd if trace is long enough
 dx=x(1+lag:lx)-x(1:lx-lag);      %displacement vector x
 dy=y(1+lag:lx)-y(1:lx-lag);      %displacement vector y
 msd=nanmean(dx.^2+dy.^2);
else
    msd=NaN;
end