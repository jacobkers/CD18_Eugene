function [firstpeaks,betterpeaks, betterpeaksvals]=JKD1_PRF_get1Dpeaksflatbottom(data,sigs,refine,plotit);
% 'Use this section for a Quicksheet'
 %------------------------------------------------------------------
    % Peakfinding routine: find peaks in a profile 'data'. 
    %A peak is a local maximum that is higher     %than a certain treshold.  
    %This treshold is found assuming a 'flat' bottom, i.e. peaks sticking
    %out of a reasonably flat base level. This base level is used to find
    %a proper measure for the standard deviation 

    %input: 1D data containing enough points to perform statistics; 
    %output: indices of the accepted peaks
    %Jacob Kers 2013
 % 'End of Quicksheet section'


if nargin<3  %DEMO
close all
data=JK00_demodata_peaks;

plotit=1;
sigs=1;
end


md=median(data);
sigma1=std(data);

localminidx=find(data(2:end-1)<data(1:end-2) &data(2:end-1)<data(3:end))+1;  %local maxes 1D
mindata=data(localminidx);
[flag,cleandata]=JKD1_PRF_outlier_flag(mindata,3,0.7,'positive',0);


sigma2=nanstd(cleandata);
md2=nanmedian(cleandata);

localmaxidx=find(data(2:end-1)>data(1:end-2) &data(2:end-1)>data(3:end))+1;  %local maxes 1D
sel=find(data(localmaxidx)>md2+sigs*sigma2);
firstpeaks=localmaxidx(sel);

%Next, perform subpixel localization
if refine
    [betterpeaks, betterpeaksvals]= Refine_Peaks(data,firstpeaks, plotit);
else
    betterpeaks=firstpeaks;
    betterpeaksvals=data(betterpeaks);
end   

if plotit %optional plot menu-------------------------------------------
    % close all;
    plot(data, 'o-'); hold on;
    plot(firstpeaks, data(firstpeaks), 'k*');
    stem(betterpeaks, betterpeaksvals, 'ko', 'MarkerFace', 'k'); %vertical lines at peak positions
end


function [betterpeaks, betterpeaksvals]= Refine_Peaks(data,firstpeaks, plotit)
%This function refines initial estimates of peak positions via parabolic
%fitting of sub-sections around each initial estimate.
%input: 'data': 1D array of points
%input: 'first peaks': array of indices of first estimates (should be integer)
%output:'betterpeaks': array of refined positions (sub-unit resolution)
%output:'betterpeaksvals': array of refined peak values (parabola maxima)
%Jacob Kers, 2013
%--------------------------------------------------------
%Run settings-----------------
subsectionsize=5; %uneven!  %3 is minimal but often OK
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
    end
end

if plotit %optional plot menu-------------------------------------------
    % close all;
    plot(data, 'o-'); hold on;
    plot(firstpeaks, data(firstpeaks), 'k*');
    stem(betterpeaks, betterpeaksvals, 'ko', 'MarkerFace', 'k'); %vertical lines at peak positions
    title('''JKD1 PRF Get1DpeaksFlatBottom'' DEMO')
end



