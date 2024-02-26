function A090_WF_SIMCorrelationDensityCurves_AbsoluteTime(batchrunindex)
 %JWJK_A:-------------------------------------------------------------------
%Correlation analysis of movie data via 1D density curves, absolute time
%
%Summary: We evaluate the changes in pattern of the chromatin over time, 
%either by comparing 1D dennity curves with each other, or the original
%images.
%
%Approach: For each  density curve, its correlation with a reference curve  
%is determined. This reference curve is chosen in two ways:
%1) ''random'': randomized: the first curve of the whole dataset is taken 
%as reference; irrespective of whether it belongs to the same cell or not. 
%We expect such correlation to be randomly fluctuating around zero for 
%most of the density curves
% 2) ''movie'': the first curve of the each cell movie is taken as reference;
% Here, we expect a non-zero correlation for early frames in the movies,
% since the patterns change not entirely from movie frame to movie frame,
%Both of these approaches can be done via density curves plotted against the 
%spatial distance axis or the genomic axis
    
%Lastly, correlation is done not via 1D density curves, but directly
%from image to image. We expect here that large-scale shape changes of the
%'donut' will be of more influence on the outcome.
%
%Input: data in .mat files stored in former analysis steps.
%
%Output: Data is presented as scatter plots and histograms. Saved are 
%Tabular excel files, .mat and summary plots.
%
%:JWJK_A-------------------------------------------------------------------
    
close all;
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);

ScatterFrameAxis=1;  %to represent curves
CorrelationModii=[{'random'} ;
                  {'movie'}];
MxFrTime=12; 
FrameAx=(1:MxFrTime)';

AvCurvesBP=zeros(MxFrTime,2);
AvCurvesDS=zeros(MxFrTime,2);


for modii=1:2
CorrelationModus=char(CorrelationModii{modii});
LabelHeaders=cell(0,0);
LabelHeaders=['Axis' 'Average' LabelHeaders];
CellNames=cell(0);

%First, make a list of cells via the first frame
sel=strfind(initval.Cell_Labels,'_t01');

for ii=1:length(sel);
    if ~isempty(sel{ii})
        CellNames=[CellNames initval.Cell_Labels{ii}];
    end
end
NCells=length(CellNames);

%two ways of correlation: 1) random: pick very first curve as density
%curve, compare all others. 2) pick first of every new cell

AllCorBP=zeros(MxFrTime,NCells+2);  AllCorBP(:,1)=FrameAx;
AllCorDS=zeros(MxFrTime,NCells+2);  AllCorDS(:,1)=FrameAx;
AllCorIm=zeros(MxFrTime,NCells+2);  AllCorIm(:,1)=FrameAx;

for cc=1:NCells
    disp(strcat('Correlation Analysis..', CorrelationModus,num2str(NCells-cc+1), 'cells to go'));
    CellName=char(CellNames(cc));
    CellName=CellName(1:end-4);
    LabelHeaders=[LabelHeaders {strcat('Cell',CellName)}];
    CellFrames=cell(0);
    
    sel=strfind(initval.Cell_Labels,CellName);
    for ii=1:length(sel);
        if ~isempty(sel{ii})
            CellFrames=[CellFrames initval.Cell_Labels{ii}];
        end
    end
    Frs=length(CellFrames); 
    CorBP_1D=zeros(MxFrTime,1)+NaN;
    CorDist_1D=zeros(MxFrTime,1)+NaN;
    CorDist_2D=zeros(MxFrTime,1)+NaN;
    
    
    for ff=1:Frs  
        cellno=char(CellFrames{ff});
        FrameTime=str2num(cellno(end-1:end));
        CellName=strcat('ResultsOfCell',cellno,'.mat'); 
        MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
        load(strcat(MatFilePath,CellName));
        if GeneralCellProps.Okayproduct
        % definition  of axes       
            DistaxisPerc=Aligned.Dist.NormAxis;
            BPaxisPerc=Aligned.BP.NormAxis;
            ContentVsDistance=Aligned.Dist.NormDensity;
            ContentVsBasePair=Aligned.BP.NormDensity; 
            ContentIm=chro_pic;
            
            
        switch CorrelationModus
            case 'movie', 
                PickRefCurve=(ff==1);
                PickExampleCurve=(ff==4);
            case 'random', 
                PickRefCurve=((ff==1)&(cc==1));
                PickExampleCurve=(ff==4);
        end
        
     
        if PickRefCurve  %pick first curve as reference
            RefCurveBP=ContentVsBasePair;
            RefCurveDist=ContentVsDistance;
            RefIm=ContentIm;
        end
        
        if PickExampleCurve  %pick first curve as reference
            ExampleCurveBP=ContentVsBasePair;
            ExampleCurveDist=ContentVsDistance;
            ExampleIm=ContentIm;
        end
        
        
        
            CorBP_1D(FrameTime)=CrossCorrelateCurves(RefCurveBP,ContentVsBasePair);
            CorDist_1D(FrameTime)=CrossCorrelateCurves(RefCurveDist,ContentVsDistance); 
            CorDist_2D(FrameTime)=CrossCorrelateImages(RefIm,ContentIm); 
            FrameAx(FrameTime)=FrameTime;
        end
    end 
    AllCorBP(:,cc+2)=CorBP_1D;
    AllCorDS(:,cc+2)=CorDist_1D;
    AllCorIm(:,cc+2)=CorDist_2D;
end

%2nd column is average
AvCorBP=(nanmean((AllCorBP(:,3:end)')))';
AvCorDS=(nanmean((AllCorDS(:,3:end)')))';
AvCorIm=(nanmean((AllCorIm(:,3:end)')))';

AllCorBP(:,2)=AvCorBP;
AllCorDS(:,2)=AvCorDS;
AllCorIm(:,2)=AvCorIm;

AvCurvesBP(:,modii)=AvCorBP;
AvCurvesDS(:,modii)=AvCorDS;
AvCurvesIm(:,modii)=AvCorIm;

Lax=length(FrameAx);
%figure(modii);
subplot(3,3,(modii-1)+1);
Scatax=(0.8*rand(Lax,1)-0.5)*ScatterFrameAxis;
plot(FrameAx+Scatax,AllCorBP(:,3:end),'ko-','MarkerSize',3); hold on;
plot(FrameAx,AllCorBP(:,2),'bo-', 'LineWidth',2);
title(strcat('genomic distance norm,',CorrelationModus));
xlabel('image index');
ylabel('correlation,n.u.');
subplot(3,3,(modii-1)+4);
plot(FrameAx+Scatax,AllCorDS(:,3:end),'ro-','MarkerSize',3); hold on;
plot(FrameAx,AllCorDS(:,2),'bo-','LineWidth',2);
title(strcat('spatial distance norm,',CorrelationModus));
xlabel('image index');
ylabel('correlation,n.u.');
pause(0.05);

subplot(3,3,(modii-1)+7);
plot(FrameAx+Scatax,AllCorIm(:,3:end),'ro-','MarkerSize',3); hold on;
plot(FrameAx,AllCorIm(:,2),'bo-','LineWidth',2);
title(strcat('whole image norm,',CorrelationModus));
xlabel('image index');
ylabel('correlation,n.u.');
pause(0.05);


%%Build and save Excel tables 
eraserarray=NaN*zeros(500,200);

%Corelation by basepair axis
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perBasePair.xlsx'),eraserarray,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perBasePair.xlsx'),LabelHeaders,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perBasePair.xlsx'),AllCorBP,'Sheet1','A2');
    save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perBasePair.mat'),'AllCorBP');


%Corelation by distance axis
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perDistance.xlsx'),eraserarray,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perDistance.xlsx'),LabelHeaders,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perDistance.xlsx'),AllCorDS,'Sheet1','A2');
save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perDistance.mat'),'AllCorDS');

%Corelation by image
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perImage.xlsx'),eraserarray,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perImage.xlsx'),LabelHeaders,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perImage.xlsx'),AllCorIm,'Sheet1','A2');
save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_Correlation',CorrelationModus,'_perImage.mat'),'AllCorIm');
end

%plot only averages
%figure(3);
subplot(3,3,3);
plot(FrameAx,AvCurvesBP,'o-', 'LineWidth',2);
title('genomic distance norm, cells av');
xlabel('image index');
ylabel('correlation,n.u.');
legend(CorrelationModii,'Location','NorthOutside','FontSize',8);
subplot(3,3,6);
plot(FrameAx,AvCurvesDS,'o-','LineWidth',2);
title('spatial distance norm, cells av');
xlabel('image index');
ylabel('correlation,n.u.');
legend(CorrelationModii,'Location','NorthOutside','FontSize',8);
pause(0.05);
subplot(3,3,9);
plot(FrameAx,AvCurvesIm,'o-','LineWidth',2);
title('whole image norm, cells av');
xlabel('image index');
ylabel('correlation,n.u.');
legend(CorrelationModii,'Location','NorthOutside','FontSize',8);

pause(0.05);





%% summary plot
figure;
subplot(2,3,1);
plot(RefCurveBP,'k-','LineWidth',2); hold on;
plot(ExampleCurveBP,'r-','LineWidth',2);
title('Reference curve, genomic axis');
xlabel('genomic position, r.u.');
ylabel('local density, r.u.');
legend('ref','example','FontSize',6,'Location','NorthOutside');

subplot(2,3,2);
plot(RefCurveDist,'k-','LineWidth',2); hold on;
plot(ExampleCurveDist,'r-','LineWidth',2);
title('Reference curve, spatial axis');
xlabel('spatial position, r.u.');
ylabel('local density, r.u.');
legend('ref','example','FontSize',6,'Location','NorthOutside');
%legend('ref','example');


subplot(4,3,3);
pcolor(RefIm); shading flat; colormap bone; axis equal; axis tight; axis off;
title('Reference image, for 2D correlation')
subplot(4,3,6);
pcolor(ExampleIm); shading flat; colormap bone; axis equal; axis tight; axis off;
title('Example image, for 2D correlation')


subplot(2,3,4);
plot(FrameAx,AvCurvesBP,'o-', 'LineWidth',2);
title('genomic distance norm, cells av');
xlabel('image index');
ylabel('correlation,n.u.');
%legend(CorrelationModii,'Location','NorthOutside','FontSize',8);
subplot(2,3,5);
plot(FrameAx,AvCurvesDS,'o-','LineWidth',2);
title('spatial distance norm, cells av');
xlabel('image index');
ylabel('correlation,n.u.');
%legend(CorrelationModii,'Location','NorthOutside','FontSize',8);
pause(0.05);
subplot(2,3,6);
plot(FrameAx,AvCurvesIm,'o-','LineWidth',2);
title('whole image norm, cells av');
xlabel('image index');
ylabel('correlation,n.u.');
%legend(CorrelationModii,'Location','NorthOutside','FontSize',8);
legend(CorrelationModii,'FontSize',8);
pause(0.05);

saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_AbsoluteTimeCorrelation.jpg'),'jpg');
saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_AbsoluteTimeCorrelation'));
    clc; 
    disp(datestr(clock));
    disp('Correlation analysis of SIM data via 1D density curves');
    disp('For each density curve, its correlation with a reference curve  is determined');
    disp('This reference curve is chosen in two ways:');
    disp('1) ''random'': randomized: the first curve of the whole dataset is taken as reference;');
    disp('irrespective of whether it belongs to the same cell or not');
    disp('We expect such correlation to be randomly fluctuating around zero for most of the density curves');
    disp('2) ''movie'': the first curve of the each cell movie is taken as reference;');
    disp('Here, we expect a non-zero correlation for early frames in the movies,');
    disp('since the patterns change not entirely from movie frame to movie frame,');
    disp('In addition, we can evaluate both type of correlation vs the spatial distance axis or the genomic axis');
    disp('Panels:');
    disp('top left:correlation per cell movie using genomic axis, all cells and average vs. frame number');
    disp('bottom left:same, using distance axis');
    disp('top middle:same as top left, randomized using one global reference curve');
    disp('bottom middle:same as bottom left, randomized using one global reference curve');    
    disp('top right:average correlation per movie and randomized, genomic axis');
    disp('bottom right: spatial axis');
    diary(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A090_DensitycurveCorrelation_Caption.txt'));

    diary off;


disp('done!');



    

function CC=CrossCorrelateCurves(RefCurve,Curve);
%this function find the symmetry center of an array

    RefCurve=RefCurve-mean(RefCurve);
   
    Curve=Curve-mean(Curve);
    
    CN=real(ifft(fft(RefCurve).*conj(fft(RefCurve))));
    CX=real(ifft(fft(Curve).*conj(fft(RefCurve))));
    CC=CX(1)/CN(1);
dum=1;

function CC=CrossCorrelateImages(RefIm,Im);
%this function find the symmetry center of an array
    %err: be sure that images have equal  size
    [rR,cR]=size(RefIm);
    [rI,cI]=size(Im);
    rU=min([rR,rI]);
    cU=min([cR,cI]);
    Im=Im(1:rU,1:cU);
    RefIm=RefIm(1:rU,1:cU);

    RefIm=RefIm/sum(RefIm(:));   
    Im=Im/sum(Im(:));  
    
    RefIm=RefIm-mean(RefIm(:));   
    Im=Im-mean(Im(:));  
    CN=real(ifft2(fft2(RefIm).*conj(fft2(RefIm))));
    CX=real(ifft2(fft2(Im).*conj(fft2(RefIm))));
    CC=CX(1)/CN(1);
dum=1;






    
    
    
    
