function [AllPeakProps,BuildCurve,AllClusterProps, AllClusters]=prf_decompose_peaks_and_clusters(PeaksCurve,Psf,showit,PeakMode)
%this function subtracts 1D Gaussians from a curve until a stop criterion
%is met. JacobKers 2016

%  AllPeakProps=[[AllPeakProps];...
%     [peakcount PeakVal Xpos Psf ThisSpotFraction(peakcount) CoveredFraction(peakcount) RelChange]];  

% Psf sets the width of Gaussians to peel off           
StopRelChange=0.03;     %0.01 means 1% of change in covered fraction
ChipIt=1;            %This is the fraction of the local maximum that ...                       %is used to build the gauss to be subtracted
MainPeakMinFract=0.20;   %this sets the treshold for 'main peaks'; the number of these (not the positions) 
%sets the number of clusters to be found by k-means clustering

                      
if nargin<4
    PeakMode='Periodic';  %'NonPeriodic';
    close all; 
    showit=1;
    Psf=23;                 %This sets the width of Gaussians to peel off  
    Pks=30;
    LL=500;
    xax=1:LL;
    PeaksCurve=0*xax;
    uxes=rand(1,Pks)*LL;
    sig=20;
    for ii=1:Pks
        ux=uxes(ii);
        PeaksCurve=PeaksCurve+OnePeak(xax,ux,sig);
    end
    %plot(xax,PeaksCurve);
end

stopit=0;
PeelCurve=PeaksCurve;
LL=length(PeaksCurve);
Xax=1:LL;


BuildCurve=0*PeaksCurve;
peakcount=0;
CoveredFraction=[];
AllPeakProps=[];


AllOnePeakCurves=[];
while ~stopit
    [PeakVal,Xpos]=max(PeelCurve);  
    dx=subpix_aroundzero(PeelCurve');
    Xpos=Xpos+dx;
    
    OnePeakCurve=ChipIt*PeakVal*OnePeak(Xax,Xpos,Psf);    
    PeelCurve=PeelCurve-OnePeakCurve;
    BuildCurve=BuildCurve+OnePeakCurve;
    peakcount=peakcount+1;
    AllOnePeakCurves(peakcount,:)=OnePeakCurve;
    CoveredFraction(peakcount)=sum(BuildCurve)/sum(PeaksCurve);
    ThisSpotFraction(peakcount)=sum(OnePeakCurve)/sum(PeaksCurve);
     
    if peakcount>1 %check progress
        RelChange=(CoveredFraction(peakcount)...
                  -CoveredFraction(peakcount-1))./...
                   CoveredFraction(peakcount);
               if (RelChange<StopRelChange) |(peakcount>50), stopit=1; end
    else
        RelChange=1;
    end    
    ResiduCurve=PeaksCurve-BuildCurve;    
    AllPeakProps=[[AllPeakProps];...
    [peakcount PeakVal Xpos Psf ThisSpotFraction(peakcount) CoveredFraction(peakcount) RelChange]];     
end

%we assume all signal goes in peaks: correction:
AllPeakProps(:,5)=AllPeakProps(:,5)/nansum(AllPeakProps(:,5));
AllPeakProps(:,6)=AllPeakProps(:,6)/(AllPeakProps(end,6));

%% 1) CLUSTER ANALYSIS: identify clusters among the main peaks. First,
%estimate a number of clusters by counting the 'main' peaks above a
%certain fraction of total content (say, 10%). 
%We find peaks by sequentially peeling off gaussians, which tends to keep
%such main peaks well separated. Then, we apply k-means
%clustering using this number. Once the clusters are identified, we can
%determine COM-position and content.
[vals,sortix]=sort(AllPeakProps(:,3));
SortPeakPosByX=AllPeakProps(sortix,:);
AllOnePeakCurvesSortByX=AllOnePeakCurves(sortix,:);
sel=find(SortPeakPosByX(:,5)>MainPeakMinFract);  %7% content peaks'

LS=max([length(sel) 1]);
Xpos=SortPeakPosByX(:,3);

[idx,C] = kmeans(Xpos,LS);
AllClusters=[];
AllClusterProps=[];
[~,LL]=size(AllOnePeakCurvesSortByX);
xax=1:LL;
for ii=1:LS
    sel2=find(idx==ii);
    if length(sel2)>1
        ThisCluster=sum(AllOnePeakCurvesSortByX(sel2,:));
    else
        ThisCluster=AllOnePeakCurvesSortByX(sel2,:);
    end
    AllClusters(ii,:)=ThisCluster;
    
    ClusterComPos=GetComposPeriodic(ThisCluster); 
    
    
    ClusterContent=sum(SortPeakPosByX(sel2,5));
    AllClusterProps(ii,:)=[ii ClusterComPos ClusterContent];
    %plot(AllClusters(ii,:)) ;
end
[~,idx]=sort(AllClusterProps(:,2));
AllClusterProps=AllClusterProps(idx,:);
AllClusterProps(:,1)=(1:LS)';
AllClusterProps(:,3)=AllClusterProps(:,3)/sum(AllClusterProps(:,3));
AllClusters=AllClusters(idx,:);
 


if showit|nargin<2
    figure;
    plot(PeaksCurve,'k-', 'LineWidth',2); hold on;
    plot(BuildCurve,'r-'); hold on;
    %plot(ResiduCurve,'k-'); hold on;
     %AllSpotProps=[peakcount PeakVal Xpos  Psf ThisSpotFraction(peakcount) CoveredFraction(peakcount) RelChange]];  
    [Pks,~]=size(AllPeakProps);
   
%      for ii=1:Pks
%          Ux=AllPeakProps(ii,3);
%          PeakVal=AllPeakProps(ii,2);
%          OnePeakCurve=PeakVal*OnePeak(Xax,Ux,Psf,PeakMode);
%          plot(OnePeakCurve); hold on;
%      end        
     plot(AllClusters', 'LineWidth',3); hold on;
     legend('original','composed', 'components');
     pause(0.5);
     %PeakVal=AllPeakProps(:,2);
     %PeakPos=AllPeakProps(:,3);
     %stem(PeakPos,PeakVal,'go','MarkerSize',8,'MarkerFaceColor','g');
     %stem(C,0*C+max(PeaksCurve),'ko','MarkerSize',8,'MarkerFaceColor','k');
       [~]=ginput(1);
       close(gcf);
end

%AllPeakProps=Eval_Peakprops(PeaksCurve,AllPeakProps,PeakMode);

function AllPeakProps=Eval_Peakprops(PeaksCurve,AllPeakProps,PeakMode)
%This function sums up some general outcomes of the cluster analysis
%output is summarized in an extra column, that allocates each component to a single 'main'cluster

%AllSpotProps=[peakcount PeakVal Xpos  Psf ThisSpotFraction(peakcount) CoveredFraction(peakcount) RelChange]];  
    LL=length(PeaksCurve);
    [LA,cc]=size(AllPeakProps);
    Psf=AllPeakProps(1,4);
    AllPeakProps(:,end+1)=zeros(LA,1);  %add 'cluster label' column
    
% identify (pre-sorted) main peaks
    Fraction=AllPeakProps(:,5);
    sel=find(Fraction>0.05);
    LS=length(sel);
    AllPeakProps(sel(1),end)=1;  %first cluster
    for jj=1:LS
        thisidx=sel(jj); 
        notself=(sel~=thisidx);
        notallocated=AllPeakProps(sel,end)==0;
        sel2=find(notself&notallocated);
        otheridxes=(sel(sel2));
        xi=AllPeakProps(thisidx,3);
        otherxes=AllPeakProps(otheridxes,3);
        dists=abs(otherxes-xi);
        if PeakMode=='Periodic'
            sel3=find(dists>LL);
            dists(sel3)=abs(LL-dists(sel3));
        end
        sel=find(dists)<2*Psf;    
        AllPeakProps(sel(jj),end)=jj;
    end
 
 function ClusterComPos=GetComposPeriodic(ThisCluster);
    LL=length(ThisCluster);
    [~,idx]=max(ThisCluster);
    xax=(1:LL)-idx+LL/2;  %peak in middle
    sel=find(xax>LL); xax(sel)=xax(sel)-LL;  %~shift
    sel=find(xax<1); xax(sel)=xax(sel)+LL;  %~shift
    %plot(xax,ThisCluster);      [~]=ginput(1); close(gcf);  
    ClusterComPos=sum(xax.*ThisCluster)/sum(ThisCluster)-LL/2+idx;
    %Need to make this periodic! (by centering on max, then correcting)
    

function PP= OnePeak(x,ux,s)
%This is the equation for a normalized1D gaussian peak value one
PP =exp (-(x-ux).^2./(2*s.^2));

function  x=subpix_aroundzero(prfx);
     %3-point subpixel fit
     xax=[-1:1:1]'; [~,mxi]=max(prfx);
     lpr=length(prfx);
     idxes=mxi-1:1:mxi+1;
     sel=find(idxes<1); idxes(sel)=idxes(sel)+lpr;
     sel=find(idxes>lpr); idxes(sel)=idxes(sel)-lpr;
     prfx=prfx(idxes);   %peak parabols with edge transfer     
     prms=polyfit(xax,prfx,2); x=-prms(2)/(2*prms(1));
     
