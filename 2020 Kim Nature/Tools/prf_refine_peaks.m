function [betterpeaks, betterpeaksvals]= prf_refine_peaks(data,firstpeaks, plotit)
%This function refines initial estimates of peak positions via parabolic
%fitting of sub-sections around each initial estimate.
%input: 'data': 1D array of points
%input: 'first peaks': array of indices of first estimates (should be integer)
%output:'betterpeaks': array of refined positions (sub-unit resolution)
%output:'betterpeaksvals': array of refined peak values (parabola maxima)
%Jacob Kers, 2013
%--------------------------------------------------------
%Run settings-----------------
subsectionsize=7; %uneven!  %3 is minimal but often OK
%------------------------------------------------

betterpeaks=0*firstpeaks;
betterpeaksvals=0*firstpeaks;
lf=length(firstpeaks);
ld=length(data);

for i=1:lf % for all initial peak positions
    %1) pick a subsection around one peak pos; ensure borders to stay within
    %data limits
    mid=firstpeaks(i);
    lo=mid-(subsectionsize-1)/2;  lo=max(lo,1);
    hi=mid+(subsectionsize-1)/2; hi=min(hi,ld);
    indexaxis=[lo:1:hi]';
    subsection=data(indexaxis);
    %2) perform a parabolic fit and optionally plot results----------------
    subs2=(subsection-nanmean(subsection))/range(subsection); %scale
    ax2=indexaxis-mean(indexaxis);
    prms=polyfit(ax2,subs2,2); %fit
    prms(1)=prms(1)*range(subsection);
    prms(2)=prms(2)*range(subsection);
    prms(3)=prms(3)*range(subsection)+nanmean(subsection); %rescale
    peakpos=-prms(2)/(2*prms(1));                       %parabol max position
    peakval=prms(1)*peakpos^2+prms(2)*peakpos+prms(3);  %parabola max value
    betterpeaks(i)=peakpos+mean(indexaxis);
    betterpeaksvals(i)=peakval;
    %3 optional plot menu-------------------------------------------------
    if plotit  
    parabolfit=prms(1)*indexaxis.^2+prms(2)*indexaxis+prms(3);
    plot(indexaxis, parabolfit, '-r');
    [~]=ginput(1);
    end
end