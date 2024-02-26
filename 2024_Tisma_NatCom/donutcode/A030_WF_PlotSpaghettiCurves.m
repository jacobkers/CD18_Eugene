function A030_WF_PlotSpaghettiCurves(batchrunindex)
%JWJK_B:-------------------------------------------------------------------
%Plotting of 1D density curves
%
%Summary: This analysis produces heat maps of the collected aligned 
%chromosome maps built by 'A020'and plots the result. further, the curves
%are grouped in 2D 'kymographs' (position-index heat maps). 
%
%Approach: The kymographs have the option to be sorted for various features. 
%Furthermore, curves can be periodically padded to highlight features near
%ori (set in the 'initval' file).
%
%Input: pre-stored .mat data.
%Output: summary plots and excel tables containing the curves.
%
%:JWJK_B-------------------------------------------------------------------

initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
BuildandSaveSpaghettiCurves('local-content',initval);
%BuildandSaveSpaghettiCurves('local-width',initval);
%BuildandSaveSpaghettiCurves('local-peakvalue',initval);


function BuildandSaveSpaghettiCurves(WhatKindOfSpaghetti,initval);

close all; pause(0.1);
initval.Cell_Labels;
allframes=length(initval.Cell_Labels);

goodcount=0;
LabelHeaders=cell(0,0);


%loop 2: paste spaghetti curves by sorting index
for jj=1:allframes
    disp(strcat('Building spaghetti curves..',num2str(allframes-jj+1),'cells to go:'));
    cellno=char(initval.Cell_Labels{jj});
    CellName=strcat('ResultsOfCell',cellno,'.mat'); 
    MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
    load(strcat(MatFilePath,CellName));               
   %% Initialize and fill two collections of chromosome density curves:   
    if jj==1     
        laa=100;
        SpaghettiVsBasePair=zeros(allframes,laa);
        SpaghettiVsDistance=zeros(allframes,laa);
        MarkerPosses=zeros(allframes,4);
    end   
    if (GeneralCellProps.Okayproduct|initval.PassAllCells)&Aligned.BP.CorrectionOK
        goodcount=goodcount+1;
        switch WhatKindOfSpaghetti
            case 'local-content'
                SpaghettiVsBasePair(goodcount,:)=Aligned.BP.NormDensity;  
                SpaghettiVsDistance(goodcount,:)=Aligned.Dist.NormDensity;
                %diagnose plot:
                if 0   
                    close all;
                   xCom=Chromosome.xCOM;
                   yCom=Chromosome.yCOM;
                   NoGaps=Chromosome.ValidIdx;
                   CX=Chromosome.CartesianContourEdge_X;
                   CY=Chromosome.CartesianContourEdge_Y;
                   MX=Chromosome.CartesianContourMax_X;
                   MY=Chromosome.CartesianContourMax_Y;
                   [ICX,ICY,OCX,OCY]=build_extra_contour_rings(CX,CY,MX,MY,xCom,yCom);
                   
                   
                   pcolor(Chromosome.picture); colormap bone; shading flat; hold on;
                    
                    subplot(2,2,1); 
                        pcolor(c2_pic); shading flat; axis equal; axis tight; hold on;
                        plot(c3.spotX, c3.spotY, 'ro','MarkerSize', 6); ; 
                        plot(xCom,yCom,'yx-', 'MarkerSize',10); hold on;
                        title('c2')
                        plot(MX(NoGaps),MY(NoGaps),'y-','LineWidth',1);
                        %plot(CX(NoGaps),CY(NoGaps),'b-','LineWidth',1);
                        plot(ICX(NoGaps),ICY(NoGaps),'b-','LineWidth',1);
                        plot(OCX(NoGaps),OCY(NoGaps),'b-','LineWidth',1);
                        plot(MX(1),MY(1),'yo'); hold off;
                    subplot(2,2,2); 
                        pcolor(c3_pic); shading flat; axis equal; axis tight; hold on;
                        plot(c3.spotX, c3.spotY, 'ro','MarkerSize', 6); hold off; 
                        title('c3')
                    subplot(2,1,2);
                    density=Aligned.Dist.NormDensity
                    true_dist_ax=Chromosome.TotalMaxPeakLength*Aligned.Dist.NormAxis/100;
                    plot(true_dist_ax,density, 'LineWidth', 2);
                    
                    title('density contour')
                    xlabel('true distance, pixels');
                    ylabel('density, a.u.');
                    [~]=ginput(1);
                    close(gcf);
                end
                
                
            case 'local-width'   
                SpaghettiVsBasePair(goodcount,:)=Aligned.BP.GeneralPropertyArray(2,:);  
                SpaghettiVsDistance(goodcount,:)=Aligned.Dist.GeneralPropertyArray(2,:);
            case 'local-peakvalue'   
                SpaghettiVsBasePair(goodcount,:)=Aligned.BP.GeneralPropertyArray(4,:);  
                SpaghettiVsDistance(goodcount,:)=Aligned.Dist.GeneralPropertyArray(4,:);
                
        end
        MarkerPosses(goodcount,:)=...
            [goodcount...
             Aligned.Orig.MarkerLabel.MarkerDistPerc...  original distance, %
             Aligned.Orig.MarkerLabel.MarkerBPPerc...    original cum.content, %
             Aligned.BP.MarkerLabel.MarkerBPPerc];       %after alignment,%
        LabelHeaders=[LabelHeaders {strcat('Cell',cellno)}];
    end
   
end
 dum=1;

%% remove empty (rejected) cells
okcells=(sum(SpaghettiVsBasePair'))>0;
SpaghettiVsBasePair=...
    SpaghettiVsBasePair(okcells,:);
SpaghettiVsDistance=...
    SpaghettiVsDistance(okcells,:);
sel=find(okcells);
LabelHeaders=LabelHeaders(sel);
LabelHeaders=['Axis' 'Average' LabelHeaders];

%% Some post-processing and definition  of axes
 DistaxisPerc=Aligned.Dist.NormAxis;
 BPaxisPerc=Aligned.BP.NormAxis;
 OriBPAxis=Aligned.Orig.NormCumDensity;
 OriDist=Aligned.Orig.NormDist;


%Pad the data front and back 25% (optional)
if initval.Padcurves>0  
    [DistaxisPerc,SpaghettiVsDistance]=...
         F008_PadData(DistaxisPerc,SpaghettiVsDistance,initval,'Xrepeat');
    [BPaxisPerc,SpaghettiVsBasePair]=...
         F008_PadData(BPaxisPerc,SpaghettiVsBasePair,initval,'Xrepeat');
end

AvSpaghettiBasepair=nanmean(SpaghettiVsBasePair);
AvSpaghettiDist=nanmean(SpaghettiVsDistance);

    
%%Build and save Excel tables (before sorting)
eraserarray=NaN*zeros(500,200);

%Using BasePair axis
SummaryTable=[BPaxisPerc' AvSpaghettiBasepair' SpaghettiVsBasePair'];
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti',WhatKindOfSpaghetti,'ByBasePair.xlsx'),eraserarray,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti',WhatKindOfSpaghetti,'ByBasePair.xlsx'),LabelHeaders,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti',WhatKindOfSpaghetti,'ByBasePair.xlsx'),SummaryTable,'Sheet1','A2');
    save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti',WhatKindOfSpaghetti,'ByBasePair.mat'),'SummaryTable');

%Using BasePair axis
SummaryTable=[DistaxisPerc' AvSpaghettiDist' SpaghettiVsDistance'];
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti',WhatKindOfSpaghetti,'ByDistance.xlsx'),eraserarray,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti',WhatKindOfSpaghetti,'ByDistance.xlsx'),LabelHeaders,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti',WhatKindOfSpaghetti,'ByDistance.xlsx'),SummaryTable,'Sheet1','A2');
save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti',WhatKindOfSpaghetti,'ByDistance.mat'),'SummaryTable');


%% sort the chromosome maps
[SpaghettiVsDistanceSort,OutSortIdx]=B050_ProcessFurther_Sort_CellMaps(SpaghettiVsDistance, initval.spaghettisortstyle);   
[SpaghettiVsBasePairSort]=B050_ProcessFurther_Sort_CellMaps(SpaghettiVsBasePair, 'Imposed', OutSortIdx);

[rr,cc]=size(SpaghettiVsBasePair);


%% Saving results
    AllCellsResultName=strcat('AllCells_',initval.expi,'.mat');
    disp('Saving collected cell results to'); disp(AllCellsResultName);
    save(strcat(initval.resultpath,AllCellsResultName),...
        'DistaxisPerc', ...
        'SpaghettiVsDistance',...
        'AvSpaghettiDist',...
        'BPaxisPerc',...
        'SpaghettiVsBasePair',...
        'AvSpaghettiBasepair');

%Plot Menu 1---------------------------------------------------------
        %% overview:
        if 0 
            figure(1);
            subplot(1,4,3); 
                pcolor(SpaghettiVsBasePairSort); colormap hot; shading flat;hold on;
                plot(MarkerPosses(:,4),MarkerPosses(:,1), 'wo');
                title('Genomic axis');
                ylabel('Cell#'); 

            subplot(1,4,4); 
                pcolor(SpaghettiVsDistanceSort); colormap hot; shading flat; hold on;
                plot(MarkerPosses(:,2),MarkerPosses(:,1), 'wo');
                title('Spatial axis');

            subplot(2,2,1);
                hold off
                [rr,cc]=size(SpaghettiVsBasePair);
                selectcurves=ceil(rr*rand(5,1));
                plot(BPaxisPerc,SpaghettiVsBasePair,'-', 'Color', [0.7 0.7 0.7]); hold on;
                plot(BPaxisPerc,SpaghettiVsBasePair(selectcurves,:),'-','LineWidth',3); hold on;
                plot(BPaxisPerc,AvSpaghettiBasepair,'ko-','MarkerSize', 4,'LineWidth',3); hold on;
                title(strcat(initval.expi,WhatKindOfSpaghetti,' vs genome %,',initval.alignmodus));
                ylabel(WhatKindOfSpaghetti); 
                xlabel('Genomic Distance, perc');
            subplot(2,2,3);
                plot(DistaxisPerc,SpaghettiVsDistance,'-', 'Color', [0.7 0.7 0.7]); hold on;
                plot(DistaxisPerc,SpaghettiVsDistance(selectcurves,:),'-','LineWidth',2); hold on;
                plot(DistaxisPerc,AvSpaghettiDist,'ko-','MarkerSize', 4,'LineWidth',2); hold off;
                title(strcat(initval.expi,WhatKindOfSpaghetti,' vs distance %,',initval.alignmodus));
                ylabel(WhatKindOfSpaghetti); 
                xlabel('Spatial distance, perc');
              saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti_',WhatKindOfSpaghetti,'curves_Overview.jpg'),'jpg');       
              plot2svg(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti_',WhatKindOfSpaghetti,'curves_Overview.svg'), gcf);
            end
        
        %% genomic axis: 
        if 0  
        figure(2); 
        subplot(1,2,1);
            plot(BPaxisPerc,SpaghettiVsBasePair,'b-','LineWidth',1); hold on;
            plot(BPaxisPerc,AvSpaghettiBasepair,'k-','MarkerSize', 4,'LineWidth',3); hold on;
            title(strcat(initval.expi,WhatKindOfSpaghetti,' vs genome percentage,',initval.alignmodus));
            ylabel(WhatKindOfSpaghetti); 
            xlabel('Genomic distance, perc');
        
        subplot(1,2,2); 
            pcolor(SpaghettiVsBasePairSort); colormap hot; shading flat;hold on;
            title('Genomic axis');
            ylabel('Cell#');    
        
        saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti_',WhatKindOfSpaghetti,'curves_Genomic.jpg'),'jpg');
        if initval.save_svg, plot2svg(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti_',WhatKindOfSpaghetti,'curves_Genomic.svg'), gcf);end
        end
        
        %% spatial axis
        if 1
            figure(3); 
            subplot(1,2,1);
                plot(DistaxisPerc,SpaghettiVsDistance,'b-','LineWidth',1); hold on;
                plot(DistaxisPerc,AvSpaghettiDist,'k-','MarkerSize', 4,'LineWidth',3); hold on;
                title(Replace_underscores(strcat(initval.expi,':',WhatKindOfSpaghetti,' vs distance %,')));
                ylabel(WhatKindOfSpaghetti); 
                xlabel('spatial distance, perc');

            subplot(1,2,2); 
                pcolor(SpaghettiVsDistanceSort); colormap hot; shading flat;hold on;
                title('spatial axis');
                ylabel('Cell#');    

            saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti_',WhatKindOfSpaghetti,'curves_Distance.jpg'),'jpg');
                if initval.save_svg
                    plot2svg(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A030_Spaghetti_',WhatKindOfSpaghetti,'curves_Distance.svg'), gcf);
                end
            end
  
        disp(['done with',WhatKindOfSpaghetti]);
        
    function ProduceFig1Aplots(DistaxisPerc,AvContentDist,OriBPAxis,OriDist)       
        figure(3);
        plot(DistaxisPerc,AvContentDist,'k-','MarkerSize', 4,'LineWidth',3); hold off;
        title('Chromosome Density vs. spatial distance from rfp');
        ylabel('Basepair density'); 
        xlabel('Spatial distance, perc');
      
        
        figure(4);
        plot(OriDist,OriBPAxis,'k-','MarkerSize', 4,'LineWidth',3); hold off;
        title('Integrated Chromosome density vs. spatial distance');
        ylabel('cumulative content, perc'); 
        xlabel('Spatial distance, perc');
        
  function [ICX,ICY,OCX,OCY]=build_extra_contour_rings(CX,CY,MX,MY,xCom,yCom);
%first, make a vector
oridist=((MY-yCom).^2+(MX-xCom).^2).^0.5;  %distance
uVX=(MX-xCom)./oridist;
uVY=(MY-yCom)./oridist;

innerdist=0.5*oridist;
outerdist=0.8*oridist;

ICX=MX-uVX.*innerdist;
ICY=MY-uVY.*innerdist;
OCX=MX+uVX.*outerdist;
OCY=MY+uVY.*outerdist;

if 0
    close all
    plot(MX,MY, 'r-'); hold on;
    %plot(CX,CY, 'b-'); 
    plot(ICX,ICY, 'm-'); 
    plot(OCX,OCY, 'm-'); hold off;
    [~]=ginput(1);
    close(gcf);
end
        

       
    


    
    
    
    
