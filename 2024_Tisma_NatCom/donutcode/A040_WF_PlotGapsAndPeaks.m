function A040_WF_PlotGapsAndPeaks(batchrunindex)
%JWJK_A:-------------------------------------------------------------------
% 
%Peak and gap analysis of density curves
%
%Summary: This function evaluates extrema of density curves. These are 
%representative for the dense clusters and open gaps between them along the
%semicircle of the chromatin pattern.
%
%Approach: minima are detected per 1D chromatin density curves. They are 
%classified as 'weak' when local density is between  50% and 100% of the
%average density (minima above average are ignored). 'Strong' minima are
%between 0 and 50%. Similarly, maxima are detected and classified. In
%addition, we evaluate how dilute or dense local chromatin is, by
%asking how much distance (in spatial units) around a minimum we should travel 
%along the chromatin semicircle for covering 2% genomic content. In the
%absence of any corrugation, this would have to be 2%. 
%
%Input: data in .mat files stored in former analysis steps.
%
%Output: Data is presented as scatter plots and histograms. Saved are 
%Tabular excel files, .mat and summary plots.
%
%:JWJK_A-------------------------------------------------------------------


close all;
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
allframes=length(initval.Cell_Labels);

AllMaxPos=[];
AllGapPos=[];
AllTerPos=[];
AllMaxVal=[];
AllGapVal=[];
AllTerVal=[];
plotGaps=[];
plotMax=[];
plotTer=[];

AllCellsGapList=[];
AllCellsPeakList=[];


for jj=1:allframes
    disp(strcat('1D Density Curve Analysis..', num2str(allframes-jj+1), 'cells to go'));
    cellno=char(initval.Cell_Labels{jj});
    NumCellLabel=BuildNumericCellLabel(cellno);
    MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
    CellName=strcat('ResultsOfCell',cellno,'.mat'); 
    load(strcat(MatFilePath,CellName)); 
    if GeneralCellProps.Okayproduct|initval.PassAllCells   
        %vs. equi-angle
        %OriBPAxis=Aligned.Orig.NormCumDensity;
        %OriDist=Aligned.Orig.NormDist;
        
        %vs. equi-distance
        Dist_Axis=Aligned.Dist.NormAxis;
        Density_Dist =Aligned.Dist.NormDensity;
        %BP_Dist=Aligned.Dist.NormCumDensity;
        BPCor_Dist=Aligned.Dist.NormCumDensityMarkercorrected;
        
        %vs. equi-content
        %BPaxisPerc=Aligned.BP.NormAxis;
        %Density_BP_=Aligned.BP.Density;

        MarkerDist=Aligned.Orig.MarkerLabel.MarkerDistPerc;  %original distance, %
        MarkerBPcor=Aligned.BP.MarkerLabel.MarkerBPPerc;       %after alignment,% 
        idx=max([round(MarkerBPcor) 1]);  %where to find it on BP percaxis
        MarkerLocalDensity=Aligned.BP.NormDensity(idx);
             
   
        %Pad the data front and back 25% (optional)
        if initval.Padcurves>0  
            [Dist_Axis,Density_Dist]=F008_PadData(Dist_Axis,Density_Dist,initval,'Xrepeat');
            [~,BPCor_Dist]=F008_PadData(Dist_Axis,BPCor_Dist,initval,'Xrepeat');
        end
    
        minima=Eval_Minima(Dist_Axis,Density_Dist,BPCor_Dist);
        minima.DistancePerc=Eval_extrema_strength(minima,Dist_Axis,Density_Dist,2);
        maxima=Eval_Maxima(Dist_Axis,Density_Dist,BPCor_Dist);
        maxima.DistancePerc=Eval_extrema_strength(maxima,Dist_Axis,Density_Dist,2);
        
        LG=length(minima.Idx);          
        thiscelllabel=zeros(LG,1)+NumCellLabel;
        
        ThisCellGaplist=...
        [thiscelllabel...
         minima.Pos'...
         minima.BPpos'...
         minima.Vals'...
         minima.DistancePerc'...
         minima.Strength'];  
          
        LG=length(maxima.Idx); 
        thiscelllabel=zeros(LG,1)+NumCellLabel;
     
        ThisCellPeaklist=...
        [thiscelllabel...
         maxima.Pos...
         maxima.BPpos'...
         maxima.Vals'...
         maxima.DistancePerc'...
         maxima.Strength'];  

        %% collect the stuff    
        plotGaps=[plotGaps;[minima.BPpos' minima.Vals']];
        plotMax=[plotMax;[maxima.BPpos',maxima.Vals']];
        plotTer=[plotTer;[MarkerBPcor MarkerLocalDensity]];  
   
        AllMaxPos=[AllMaxPos maxima.BPpos];
        AllGapPos=[AllGapPos minima.BPpos];
        AllTerPos=[AllTerPos MarkerBPcor];
        
        %AllMaxVal=[AllMaxVal maxima.Vals];
        %AllGapVal=[AllGapVal minima.Vals];
        AllMaxVal=[AllMaxVal maxima.DistancePerc];
        AllGapVal=[AllGapVal minima.DistancePerc];        
        
        AllTerVal=[AllTerVal MarkerLocalDensity];
        
        %add the new gap and peak entries
        AllCellsGapList=[AllCellsGapList; ThisCellGaplist];
        AllCellsPeakList=[AllCellsPeakList; ThisCellPeaklist];
        
        %add the ter info
        TerEntry=[NumCellLabel MarkerDist MarkerBPcor MarkerLocalDensity NaN NaN];
        AllCellsGapList=[AllCellsGapList; TerEntry];
        AllCellsPeakList=[AllCellsPeakList; TerEntry];        
    end
end

%% scatter plot gaps and peaks
    if 1
        subplot(2,2,1);
        %figure(1);

        plot(plotGaps(:,1),plotGaps(:,2),'ko','MarkerSize',3); hold on;
        plot(plotMax(:,1),plotMax(:,2),'ko','MarkerSize',3,'MarkerFaceColor','y'); hold on;
        plot(plotTer(:,1),plotTer(:,2), 'bo','MarkerSize',3, 'MarkerFaceColor','b'); hold on;
        title(strcat(initval.expi,':location of maxima and minima'));

        legend('gaps','peaks','marker')
        xlabel('genomic position, %');
        ylabel('local density, %');
    end
    
    %Save the summaries
    ColNames=[, {'label'}, {'position, distance'} ,{'position, basepair'}, {'value'},...
           {'distance for 2%genome'}, {'strength'}];
       eraserarray=NaN*zeros(500,200); 
        xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensityCurvesGapData.xlsx'),eraserarray,'Sheet1','A1');
        xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensityCurvesGapData.xlsx'),ColNames,'Sheet1','A1');
        xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensityCurvesGapData.xlsx'),AllCellsGapList,'Sheet1','A2');
            save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensityCurvesGapData.mat'),'AllCellsGapList');
        xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensityCurvesPeakData.xlsx'),eraserarray,'Sheet1','A1');
        xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensityCurvesPeakData.xlsx'),ColNames,'Sheet1','A1');
        xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensityCurvesPeakData.xlsx'),AllCellsPeakList,'Sheet1','A2');
            save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensityCurvesPeakData.mat'),'AllCellsPeakList');
     
    
            
            
%% histogram of peaks and gaps vs BP positions
    %figure(2);
    subplot(2,2,2);
    Build_PositionHistograms(AllMaxPos,AllGapPos,AllTerPos,initval);
  
%%  %% build a histogram of values around ori and ter
    %figure(3);
    subplot(2,2,3);
    Build_NearLabelPositionHistograms(AllMaxPos,AllGapPos,AllTerPos,initval);
   

%% classify the strengt of gaps (BP distance they cover 
    %figure(4);
    subplot(2,2,4);
    Build_ValueHistograms(AllGapPos,AllGapVal,40,initval);
    
%% output summary
    saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensitycurveAnalysis.jpg'),'jpg');
    saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensitycurveAnalysis'));
    clc; 
    diary(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A040_DensitycurveAnalysis_Caption.txt'));
    disp(datestr(clock));
    disp('Analysis of peaks and gaps via 1D density curves');
    disp('top left:scatter plot of local density at peaks, gaps and cfp labels; 1% is average density');
    disp('top right:same, counts per position');
    disp('bottom left:occurence vs. distance a gap of peak is found from the cfp position');
    disp('bottom right:strength of gaps (expressed as distance around gaps that covers 2% of genomic content');   
    diary off;
    disp('done');
 
   
function Build_ValueHistograms(AllGapPos,AllGapVal,span,initval);

    selA=find((AllGapPos>83.3)|(AllGapPos<16.7));  %33.3%near ori
    selB=find((AllGapPos<66.7)&(AllGapPos>33.3));  %33%near ter
    selC=find((AllGapPos>16.7)&(AllGapPos<33.3)|...
              (AllGapPos>66.7)&(AllGapPos<83.3));  %33%rest
    AllGapValA=AllGapVal(selA);
    AllGapValB=AllGapVal(selB);
    AllGapValC=AllGapVal(selC);

    binax=linspace(0, span,20);
    GapValHistA=hist(AllGapValA,binax);
    GapValHistB=hist(AllGapValB,binax);
    GapValHistC=hist(AllGapValC,binax);
    shft=0.3;
    bar(binax+shft,GapValHistB,'b','BarWidth',0.5); hold on;
    bar(binax,GapValHistA,'r','BarWidth',0.5); hold on;
    bar(binax-shft,GapValHistC,'w','BarWidth',0.5); hold on;
    legend('Gaps near cfp','Gaps near rfp', 'Other Gaps');
    xlabel('distance covered,%');
    ylabel('counts');


 function Build_PositionHistograms(AllMaxPos,AllGapPos,AllTerPos,initval)
    skp=4;
    binax=-initval.Padcurves:skp:100+initval.Padcurves;
    MaxHist=hist(AllMaxPos,binax);
    GapHist=hist(AllGapPos,binax); 
    TerHist=hist(AllTerPos,binax);
    MaxHist(end)=0;
    GapHist(end)=0;
    TerHist(end)=0;
    
    shft=1/3;
    bar(binax-skp*shft,MaxHist,'y', 'BarWidth',1.1*shft); hold on;
    bar(binax,GapHist,'k','BarWidth',1.1*shft); hold on;
    bar(binax+skp*shft,TerHist,'b','BarWidth',1.1*shft); hold on;
    title('gap and peak local occurence');
    legend('peaks','gaps','marker');
    xlabel('genomic position, %');
    ylabel('counts');
    
 function Build_NearLabelPositionHistograms(AllMaxPos,AllGapPos,AllTerPos,initval)
    skp=1;     binax=1:3:25;
    
    %first, set absolute position
    AllMaxPos=abs(AllMaxPos);
    AllGapPos=abs(AllGapPos);
    
    MaxHist=hist(AllMaxPos,binax);
    GapHist=hist(AllGapPos,binax); 

    MaxHist(end)=0;
    GapHist(end)=0;

    shft=1/2;
    bar(binax-skp*shft,MaxHist,'y', 'BarWidth',1.1*shft); hold on;
    bar(binax,GapHist,'k','BarWidth',1.1*shft); hold on;
    title('gap and peak relative to marker');
    legend('peaks','gaps')
    xlabel('genomic distance from marker, %');
    ylabel('counts');
    

    
 
    



   
