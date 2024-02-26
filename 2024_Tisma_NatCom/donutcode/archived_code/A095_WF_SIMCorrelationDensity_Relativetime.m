function A095_WF_SIMCorrelationDensity_Relativetime(batchrunindex)
 %JWJK_A:-------------------------------------------------------------------
%Correlation analysis of movie data via 1D density curves, relative time
%
%Summary: We evaluate the changes in pattern of the chromatin, 
%either by comparing 1D density curves with each other, or the original
%images. Changes are evalauted as a function of relative time, i.e. the
%frame time difference between the two curves or images.
%
%Approach: For each frame time difference, corresponding density curves are correlated.
%The results are grouped in two ways: 
% 1) ''diff cell'': randomized: density cureves of different cells
% We expect such correlation to be randomly fluctuating around zero for 
% most of the density curves
% 2) ''same cell'': 
% Here, we expect a non-zero correlation samll frame diffferences,
% since the patterns change not entirely from movie frame to movie frame,
% In addition, we can evaluate both type of correlation vs the 
% spatial distance axis or the genomic axis
    
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
limitedset=1E6;
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);

ScatterFrameAxis=1;  %to represent curves
CorrelationModii=[{'random'} ;
                  {'movie'}];

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
NCells=min([length(CellNames) limitedset]);

CorrelationResults=zeros(NCells*NCells,10);  %we'll crop later
entrycounter=0;
%format output: index mv1 fr1 mv2 fr2 cor1D_BP cor1D_DS cor2D samecell 

%define example curves
exmv1=1; 
exmv2=2;
exfr1=1; 
exfr2=2;


for cc1=1:NCells  %for all cells:        
    disp(strcat('Correlation Analysis..',num2str(NCells-cc1+1), 'cells to go'));
    CellName=char(CellNames(cc1));
    CellName=CellName(1:end-4);
    CellFrames=cell(0);    
    sel=strfind(initval.Cell_Labels,CellName);
    for ii=1:length(sel);
        if ~isempty(sel{ii})
            CellFrames=[CellFrames initval.Cell_Labels{ii}];
        end
    end
    Frs1=length(CellFrames);   
    for ff1=1:Frs1  %for all frames of a cell
        cellno1=char(CellFrames{ff1});
        FrameTime=str2num(cellno1(end-1:end));
        CellName1=strcat('ResultsOfCell',cellno1,'.mat'); 
        MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
        load(strcat(MatFilePath,CellName1));
        if GeneralCellProps.Okayproduct
            %load the info of the 'primary frame'
            DistaxisPerc1=Aligned.Dist.NormAxis;
            BPaxisPerc1=Aligned.BP.NormAxis;
            ContentVsDistance1=Aligned.Dist.NormDensity;
            ContentVsBasePair1=Aligned.BP.NormDensity; 
            ContentIm1=chro_pic;            
            %load all cells of all movies (includes self) higher time frames
            %remember some example curves for later          
            for cc2=1:NCells  %again, for all cells
                CellName2=char(CellNames(cc2));
                CellName2=CellName2(1:end-4);
                CellFrames2=cell(0);    
                sel2=strfind(initval.Cell_Labels,CellName2);
                for ii=1:length(sel2);
                    if ~isempty(sel2{ii})
                        CellFrames2=[CellFrames2 initval.Cell_Labels{ii}];
                    end
                end
                Frs2=length(CellFrames2);   
                for ff2=1:Frs2
                    if ff2>=ff1 %only equal or higher times
                        %load the properties of the second cell frame
                        cellno2=char(CellFrames2{ff2});
                        FrameTime2=str2num(cellno2(end-1:end));
                        CellName2=strcat('ResultsOfCell',cellno2,'.mat'); 
                        MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
                        load(strcat(MatFilePath,CellName2));
                        if GeneralCellProps.Okayproduct
                            entrycounter=entrycounter+1;
                            %load the info of the 'secondary frame'
                            DistaxisPerc2=Aligned.Dist.NormAxis;
                            BPaxisPerc2=Aligned.BP.NormAxis;
                            ContentVsDistance2=Aligned.Dist.NormDensity;
                            ContentVsBasePair2=Aligned.BP.NormDensity; 
                            ContentIm2=chro_pic;
                            
                            CorBP_1D=CrossCorrelateCurves(ContentVsBasePair1,ContentVsBasePair2);
                            CorDist_1D=CrossCorrelateCurves(ContentVsDistance1,ContentVsDistance2); 
                            CorDist_2D=CrossCorrelateImages(ContentIm1,ContentIm2); 
                            samecell=1.0*(cc1==cc2);
                            deltaframe=ff2-ff1;
                            CorrelationResults(entrycounter,:)=...
                            [entrycounter cc1 ff1 cc2 ff2 CorBP_1D CorDist_1D CorDist_2D deltaframe samecell];    
                            if (cc1==exmv1 & samecell)
                                if ff1==exfr1; %example curves 1
                                    RefCurveDist=ContentVsDistance1;
                                    RefCurveBP=ContentVsBasePair1;    
                                    RefIm=ContentIm1;
                                end
                                if ff2==exfr2; %example curves 2
                                    ExampleCurveDist=ContentVsDistance2;
                                    ExampleCurveBP=ContentVsBasePair2;
                                    ExampleIm=ContentIm2;
                                end
                            end
                        end
                    end
                end 
            end            
        end
    end 
end
CorrelationResults=CorrelationResults(1:entrycounter,:);  %crop unused rows

%% condensing the results per delta-time
%[entrycounter cc1 ff1 cc2 ff2 CorBP_1D CorDist_1D CorDist_2D deltaframe samecell];
MaxDifT=max(CorrelationResults(:,9));
CorrelationPerDT=zeros(MaxDifT+1,25);
%deltaframe 
for dfr=0:MaxDifT    
    CorBP_1D_sameCell=CondenseCorrelationResults(CorrelationResults,dfr,6,1);
    CorDist_1D_sameCell=CondenseCorrelationResults(CorrelationResults,dfr,7,1); 
    CorDist_2D_sameCell=CondenseCorrelationResults(CorrelationResults,dfr,8,1); 
    
    CorBP_1D_diffCell=CondenseCorrelationResults(CorrelationResults,dfr,6,0);
    CorDist_1D_diffCell=CondenseCorrelationResults(CorrelationResults,dfr,7,0); 
    CorDist_2D_diffCell=CondenseCorrelationResults(CorrelationResults,dfr,8,0); 
    
    CorrelationPerDT(dfr+1,:)=[dfr,...
    CorBP_1D_sameCell, CorDist_1D_sameCell, CorDist_2D_sameCell,...
    CorBP_1D_diffCell, CorDist_1D_diffCell, CorDist_2D_diffCell];
end

    ColNames=[
        {'delta-frame'} , ...
        {'number found'}, ...
        {'CorBP_1D_sameCell_AV'}, ...
        {'CorBP_1D_sameCell_SEM'}, ...
        {'CorBP_1D_sameCell_STD'},...
        {'number found'}, ...
        {'CorDS_1D_sameCell_AV'}, ...
        {'CorDS_1D_sameCell_SEM'}, ...
        {'CorDS_1D_sameCell_STD'},...
        {'number found'}, ...
        {'Cor2D_sameCell_AV'}, ...
        {'Cor2D_sameCell_SEM'}, ...
        {'Cor2D_sameCell_STD'},...
        {'number found'}, ...
        {'CorBP_1D_diffCell_AV'}, ...
        {'CorBP_1D_diffCell_SEM'}, ...
        {'CorBP_1D_diffCell_STD'},...
        {'number found'}, ...
        {'CorDS_1D_diffCell_AV'}, ...
        {'CorDS_1D_diffCell_SEM'}, ...
        {'CorDS_1D_diffCell_STD'},...
        {'number found'}, ...
        {'Cor2D_diffCell_AV'}, ...
        {'Cor2D_diffCell_SEM'}, ...
        {'Cor2D_diffCell_STD'},...
        ];   
    
    %%Build and save Excel tables 
    eraserarray=NaN*zeros(NCells*NCells,30);
    
    %Corelation by basepair axis
    xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A095_CorrelationDeltaFrame.xlsx'),eraserarray,'Sheet1','A1');
    xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A095_CorrelationDeltaFrame.xlsx'),ColNames,'Sheet1','A1');
    xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A095_CorrelationDeltaFrame.xlsx'),CorrelationPerDT,'Sheet1','A2');
    save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A095_CorrelationDeltaFrame.mat'),'CorrelationPerDT');

    dum=1;
    if 1   
    
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

FrameAx=CorrelationPerDT(:,1);
AvCurvesBP_same=CorrelationPerDT(:,3);
AvCurvesDS_same=CorrelationPerDT(:,7);
AvCurvesIm_same=CorrelationPerDT(:,11);
AvCurvesBP_diff=CorrelationPerDT(:,15);
AvCurvesDS_diff=CorrelationPerDT(:,19);
AvCurvesIm_diff=CorrelationPerDT(:,23);

subplot(2,3,4);
plot(FrameAx,AvCurvesBP_same,'bo-', 'LineWidth',2); hold on;
plot(FrameAx,AvCurvesBP_diff,'go-', 'LineWidth',2); hold on;
title('genomic distance norm, cells av');
xlabel('frame difference');
ylabel('correlation,n.u.');

subplot(2,3,5);
plot(FrameAx,AvCurvesDS_same,'bo-','LineWidth',2); hold on ;
plot(FrameAx,AvCurvesDS_diff,'go-','LineWidth',2); hold on;
title('spatial distance norm, cells av');
xlabel('frame difference');
ylabel('correlation,n.u.');
legend('same cell', 'different','FontSize',8);

pause(0.05);
subplot(2,3,6);
plot(FrameAx,AvCurvesIm_same,'bo-','LineWidth',2); hold on;
plot(FrameAx,AvCurvesIm_diff,'go-','LineWidth',2); hold on;
title('whole image norm, cells av');
xlabel('frame difference');
ylabel('correlation,n.u.');
pause(0.05);

saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A095_RelativeCorrelation.jpg'),'jpg');
saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A095_RelativeCorrelation'));
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
end

function CorSummary=CondenseCorrelationResults(CorrelationResults,deltafr,col,samecell);
%as title says
%in: [entrycounter cc1 ff1 cc2 ff2 CorBP_1D CorDist_1D CorDist_2D deltaframe samecell];
    sel=find((CorrelationResults(:,9)==deltafr)&(CorrelationResults(:,10)==samecell));
    Nfound=length(sel);
    Cor_AV=nanmean(CorrelationResults(sel,col));
    Cor_STD=nanstd(CorrelationResults(sel,col));
    Cor_SEM=2*Cor_STD/Nfound^0.5; %2sigma
    CorSummary=[Nfound Cor_AV Cor_SEM Cor_STD];

function CC=CrossCorrelateCurves(RefCurve,Curve);
%this function find the symmetry center of an array
    RefCurve=RefCurve-mean(RefCurve);   
    Curve=Curve-mean(Curve);    
    CN=real(ifft(fft(RefCurve).*conj(fft(RefCurve))));
    CX=real(ifft(fft(Curve).*conj(fft(RefCurve))));
    CC=CX(1)/CN(1);

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






    
    
    
    
