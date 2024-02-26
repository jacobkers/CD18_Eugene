function [AllPeakProps,BuildCurve,AllClusterProps, AllClusters]=prf_decompose_peaks_and_clusters(PeaksCurve,Psf,merge_option,showit)
%JWJK_C*:-------------------------------------------------------------------
%Title: decompose a curve in peaks by subtraction
%Summary: This function subtracts 1D Gaussians from a curve until a stop criterion
%is met. JacobKers 2016
%Input: input curve, estimated point spread function, show-flag. Runs also in demo-mode 
%Output: properties and curves
%[peakcount PeakVal Xpos Psf ThisSpotFraction(peakcount) CoveredFraction(peakcount) RelChange]];
%References: Jacob Kers 2016
%:JWJK_C*-------------------------------------------------------------------

%settings:
%  AllPeakProps=[[AllPeakProps];...
%  [peakcount PeakVal Xpos Psf ThisSpotFraction(peakcount) CoveredFraction(peakcount) RelChange]];  
% Psf sets the width of Gaussians to peel off           
StopRelChange=0.1;     %0.01 means 1% of change in covered fraction
ChipIt=1;            %This is the fraction of the local maximum that ...                       %is used to build the gauss to be subtracted
MainPeakMinFract=0.08;  

                     
if nargin<3  %DEMO mode
    close all;     showit=1;
    merge_option='merge_by_distance' ; %'kmeans';
    %merge_option='kmeans';
    Psf=20;  Pks=12;  LL=500;    xax=1:LL;   
    sig=20; %used to build
    PeaksCurve=0*xax; uxes=rand(1,Pks)*LL;    
    for ii=1:Pks
        ux=uxes(ii);
        PeaksCurve=PeaksCurve+prf_one_gauss_peak(xax,ux,sig,0);
    end
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
    dx=prf_subpix_aroundzero(PeelCurve);
    Xpos=Xpos+dx;
    
    OnePeakCurve=ChipIt*PeakVal*prf_one_gauss_peak(Xax,Xpos,Psf,0);    
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

%% 1) CLUSTER ANALYSIS: identify clusters among the main peaks. 
%1) option 'kmeans' 
%First,estimate a number of clusters by counting the 'main' peaks above a
%certain fraction of total content (say, 10%). 
%We find peaks by sequentially peeling off gaussians, which tends to keep
%such main peaks well separated. Then, we apply k-means
%clustering using this number. Once the clusters are identified, we can
%determine COM-position and content.
%2) option merge_by_distance. For irregular but otherwise well-separated
%blobs, we use a simple distance criterion to collect adjacent 'chains'
switch merge_option
    case 'kmeans'
        [AllClusterProps,AllClusters]=get_clusters_by_kmeans(AllOnePeakCurves,AllPeakProps,MainPeakMinFract);
     case 'merge_by_distance'
         [AllClusterProps,AllClusters]=get_clusters_by_distance(AllOnePeakCurves,AllPeakProps,Psf);         
end
LS=length(AllClusterProps(:,1));
[~,idx]=sort(AllClusterProps(:,2));
AllClusterProps=AllClusterProps(idx,:);
AllClusterProps(:,1)=(1:LS)';
AllClusterProps(:,3)=AllClusterProps(:,3)/sum(AllClusterProps(:,3));
AllClusters=AllClusters(idx,:);
 
if showit;
    %figure;
    plot(PeaksCurve,'k-', 'LineWidth',2); hold on;
    plot(BuildCurve,'r-'); hold on;
    [Pks,~]=size(AllPeakProps);
         
     plot(AllClusters', 'LineWidth',3); hold on;
     pause(0.1);
     [~]=ginput(1);
end

 function [AllClusterProps,AllClusters]=get_clusters_by_kmeans(AllOnePeakCurves,AllPeakProps,MainPeakMinFract);
 %Build clusters by first finding brighmost peaks, use this as clusternumber, ...
 %then apply kmeans for this number. 
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
            ClusterComPos=prf_get_periodic_com(ThisCluster);   
            ClusterContent=sum(SortPeakPosByX(sel2,5));
            AllClusterProps(ii,:)=[ii ClusterComPos ClusterContent];
            %plot(AllClusters(ii,:)) ;
        end
 
function [AllClusterProps,AllClusters]=get_clusters_by_distance(AllOnePeakCurves,AllPeakProps,Psf);
        %collect spots within miimal distance of each other
         xposses=(AllPeakProps(:,3)); Lx=length(xposses);
         [xposses,idx]=sort(xposses); 
         AllPeakProps=AllPeakProps(idx,:); 
         AllOnePeakCurves=AllOnePeakCurves(idx,:);
          if 0
            %figure;
             plot(AllOnePeakCurves','g--', 'LineWidth',1); hold on;
             plot(sum(AllOnePeakCurves),'k-','LineWidth',2); hold on;
             pause(0.1);
          end    
         %get distance matrix
         dist_matrix=((abs(xposses-xposses')<2*Psf)&(abs(xposses-xposses')>0));
         allmerged=0; 
         clusterlabels=0*xposses;  clusterlabels(1)=1; 
         thiscluster=1; %first cluster:first spot
         while allmerged==0             
            no_change=0;  
            while ~no_change
               spots_of_thisclus=find(clusterlabels==thiscluster);
               leftspots=find(clusterlabels==0);  %uncategorized
               LC=length(spots_of_thisclus);
               for ispot=1:LC  %for all spots in this cluster
                   xs=xposses(spots_of_thisclus(ispot));
                   i_otherspots=find(abs(xposses(leftspots)-xs)<2.5*Psf) ;
                   if ~isempty(i_otherspots)
                        new_ones=leftspots(i_otherspots);
                        clusterlabels(new_ones)=thiscluster;
                   end
                   dum=1;
               end
               dist_matrix;
               LC_nw=length(find(clusterlabels==thiscluster));
               if LC_nw==LC  %no change, finalize and go to next free spot (if available)
                   no_change=1;
                   
                   %'cluster complete'
                   %finalize
                   i_final_spots_for_this_cluster=(find(clusterlabels==thiscluster));
                   ThisCluster=sum(AllOnePeakCurves(i_final_spots_for_this_cluster,:),1);
                   AllClusters(thiscluster,:)=ThisCluster;  
                   ClusterComPos=prf_get_periodic_com(ThisCluster);   
                   ClusterContent=sum(AllPeakProps(i_final_spots_for_this_cluster,5));
                   AllClusterProps(thiscluster,:)=[thiscluster ClusterComPos ClusterContent];
                   if 0
                       %close all
                       plot(AllClusters', 'LineWidth',2) ;                
                   end    
                   
                   %go to next cluster
                   thiscluster=thiscluster+1;
                   leftspots=find(clusterlabels==0);  %uncategorized
                   if ~isempty(leftspots)
                        clusterlabels(leftspots(1))=thiscluster;  %first free spot
                   else
                       allmerged=1;
                       no_change=1;
                   end
               end
             end
         end
         
         
         allmerged=1;
            
    


     
