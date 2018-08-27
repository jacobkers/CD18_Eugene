function [peakprops,buildcurve]=peel_peaks_from_profile(PeaksCurve,Psf,showit)
%this function subtracts 1D Gaussians from an image until a stop criterion
%is reached. JacobKers 2016

MaxPeakNo=7;           
StopRelChange=0.2;     %0.01 means 1% of change in covered fraction
ChipIt=1;            %This is the fraction of the local maximum that ...                        %is used to build the gauss to be subtracted

if nargin<2                %This sets the width of Gaussians to peel off   
    close all; 
    showit=1;
    Psf=23;     %This sets the width of Gaussians to peel off              
    Pks=20;
    LL=500;
    xax=1:LL;
    PeaksCurve=0*xax;
    uxes=rand(1,Pks)*LL;
    sig=20;
    for ii=1:Pks
        ux=uxes(ii);
        PeaksCurve=PeaksCurve+OnePeak(xax,ux,sig);
    end
    plot(xax,PeaksCurve);
end

stopit=0;
PeelCurve=PeaksCurve;
LL=length(PeaksCurve);
Xax=1:LL;
buildcurve=0*PeaksCurve;
peakcount=1;
CoveredFraction=[];
peakprops=[];

close all;
while ~stopit
    [PeakVal,Xpos]=max(PeelCurve);                 
    OnePeakCurve=ChipIt*PeakVal*OnePeak(Xax,Xpos,Psf);    
    PeelCurve=PeelCurve-OnePeakCurve;
    buildcurve=buildcurve+OnePeakCurve;
    peakcount=peakcount+1;
    if peakcount>1
    CoveredFraction(peakcount)=sum(buildcurve)/sum(PeaksCurve);
    ThisSpotFraction(peakcount)=sum(OnePeakCurve)/sum(PeaksCurve);     
    RelChange=(CoveredFraction(peakcount)...
              -CoveredFraction(peakcount-1))./...
               CoveredFraction(peakcount);
   stopit=( (RelChange<StopRelChange)| ...
            (peakcount>MaxPeakNo)&...
            (peakcount>=1));
    RelChange=(RelChange*(peakcount>1)+...;
              1.0*(peakcount<=1));
    end
    ResiduCurve=PeaksCurve-buildcurve;    
    peakprops=[[peakprops];...
    [peakcount PeakVal Xpos  Psf ThisSpotFraction(peakcount) CoveredFraction(peakcount) RelChange]];     
end


if showit|nargin<2
    plot(PeaksCurve,'k-', 'LineWidth',2); hold on;
    plot(buildcurve,'r-'); hold on;
    plot(ResiduCurve,'k-'); hold on;
     %AllSpotProps=[peakcount PeakVal Xpos  Psf ThisSpotFraction(peakcount) CoveredFraction(peakcount) RelChange]];  
    [Pks,~]=size(peakprops);
     for ii=1:Pks
         Ux=peakprops(ii,3);
         PeakVal=peakprops(ii,2);
         OnePeakCurve=PeakVal*OnePeak(Xax,Ux,Psf);
         plot(OnePeakCurve); hold on;
     end    
     legend('original','composed','residu', 'components');
     [~]=ginput(1);
end
   

function PP= OnePeak(x,ux,s)
%This is the equation for a normalized1D gaussian peak value one
PP =exp (-(x-ux).^2./(2*s.^2));




