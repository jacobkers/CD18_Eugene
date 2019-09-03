function [msdcurve,allpoints]=prf_mean_square_displacement(x,y)
% %This function calculates the mean square displacement of an x or xy-trajectory
% %afo lag time in frames. output is a  curve
%JacobKers 2019------------------------------------------------------------

%% demo mode
if nargin==1
    y=0*x;
end
if nargin<1
    close all;
    x=cumsum(randn(1000,1));
    y=cumsum(randn(1000,1));
end

%%
lx=length(x);
allpoints=[];
maxlag=lx/2;
for lag=1:maxlag-1
    if lag<lx/2                    %only calculate msd if trace is long enough
        dx=x(1+lag:lx)-x(1:lx-lag);      %displacement vector x
        dy=y(1+lag:lx)-y(1:lx-lag);      %displacement vector y
        msdcurve(lag)=nanmean(dx.^2+dy.^2);
        allpoints=[allpoints ; [lag+0*dx   dx.^2+dy.^2]];
    else
        msdcurve(lag)=NaN;
    end
end

%% demo mode
if nargin<1
    subplot(1,2,1);
        plot(x); hold on;
        plot(y);
        legend('x-trace','y-trace');
        title('trajectory')
    subplot(1,2,2);
    plot(msdcurve,'r-');
    title('Mean square displacement curve')
    ylabel('msd, L^2');
    xlabel('lag time, time units');
end