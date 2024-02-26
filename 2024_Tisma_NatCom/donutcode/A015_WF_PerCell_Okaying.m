function A015_WF_PerCell_Okaying(batchrunindex,dryrun)
%JWJK_A:-------------------------------------------------------------------
%Title: 
%Automatic multi-parameter cell screening
%Summary: This program performs automatic screening on proper cell size, 
%ori-ter positions etc. reports the results to the command line and saves 
%a summary plot and excel table.   
%Approach: the chromatin structure, label positions and cell wall positions 
%are used to evaluate to what extend the cell has a well-analyzable donut. 
%(in terms of size and center hole depth) and labels in a position that is roughly to be
%expected based on their genomic location.
%Input: .mat files from former analysis steps.
%Output: 
%1)a structure 'summary' containing various rejection parameters     
%initialcount: counter to see how much cells are accepted  
%cellsizeOK:  Cell size > 15, std <8 pixels
%chromosomesizeOK: Chromosome size > 5, std <3.5 pixels
%donutOK: center darker dan 0.3 of chromosome max, (0.6 for SIM)
%teroriangleOK: ter-ring center-ori angle large than 60
%terorilociiOK: ori and ter at least 0.3 times average
%chromosome radius away from center
        %totalok: Product of all filters
%2)Tabular data excel, .mat and summary plots
%3)A directory with only 'accepted' summary pictures of cells.
%:JWJK_A-------------------------------------------------------------------

actions.do_the_cell_jpg=1 ; %if you just want the bar plots, set to 0

initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
allframes=length(initval.Cell_Labels);
close all;

imoutdir=strcat(initval.resultpath,'CellImages_Screened',initval.DirSep);
imoutdir_rej=strcat(initval.resultpath,'CellImages_Rejected',initval.DirSep);

if 1
    if isdir(imoutdir), rmdir(imoutdir,'s'); end
    if isdir(imoutdir_rej), rmdir(imoutdir_rej,'s'); end
    mkdir(imoutdir); 
    mkdir(imoutdir_rej); 
end


%find files listing rejected cells
UserRejectedCellsList=GetRejectedCellNumbers(initval);


%counter to see how much cells are accepted
Summary.initialcount=allframes;  
Summary.cellsizeOK=0;  
Summary.chromosomesizeOK=0;
Summary.donutOK=0;
Summary.teroriangleOK=0;
Summary.terorilociiOK=0;
Summary.totalok=0;
Summary.spotok=0;
Summary.userok=0;
SummaryTable=[];
%contains: index cellnumber 
            %[cell size: mean sd  /ok] 
            %[chro size: mean std /ok]
            %[donutness: valuemax valuecenter /ok]
            %[ter-center-ori-angle: value /OK]
            %[ter-ori-locii: value: /OK]
            %allOK
            
            
%screening run: collect each property that is screened to explicitly show the cutoffs
summary_vals=struct('meancellradius',[],...
                    'stdcellradius', [],...
                    'cell_circularity', [],...
                    'averageradius_chro', [],...
                    'stdradius_chro', [],...
                    'chro_circularity', [],...
                    'donutness', [],...
                    'ori_ter_angle', [],...
                    'distance_ter_center', [],...
                    'distance_ori_center', []);
                

for jj=1:allframes  
    GeneralCellProps=struct('OkayCell','0');
    cellno=char(initval.Cell_Labels{jj});
    disp(strcat('Screening..',num2str(allframes-jj+1),'cells to go'));
    CellName=strcat('ResultsOfCell',cellno,'.mat'); 
    MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
    if ~ dryrun
    %% load the chromosome and cell info
    load(strcat(MatFilePath,CellName),'Chromosome','c4','c3', 'c4_xt', 'c4_pic','c3_pic','c4_pic_xt','c3_pic_xt');  
    
    cellno=char(initval.Cell_Labels{jj});
    NumCellLabel=BuildNumericCellLabel(cellno);
    CellSpecs=[jj NumCellLabel];
        
    
    %% Filter0: UserReject
    sel=find(UserRejectedCellsList==CellSpecs(2));
    if ~isempty(sel)
        GeneralCellProps.OkayByUser=0;
    else
        GeneralCellProps.OkayByUser=1;
    end
    Summary.userok=Summary.userok+GeneralCellProps.OkayByUser;    
    
    %% Filter1: Cell size
    celledge = Chromosome.RadialCellEdge; %252 Vector
    NAN = isnan(celledge); % Since there are some 'NAN'(not a number) within the matrix/vector
    validindex = NAN==0;
    celledgevalid = celledge(validindex);
    meancellradius = mean(celledgevalid);
    stdcellradius = std(celledgevalid);
    cell_contourlength=nansum(((Chromosome.CartesianCellwallX(2:end)-Chromosome.CartesianCellwallX(1:end-1)).^2+...
                            (Chromosome.CartesianCellwallY(2:end)-Chromosome.CartesianCellwallY(1:end-1)).^2).^0.5);  
    cell_circularity=cell_contourlength/meancellradius/(2*pi);
    %evaluate:
    GeneralCellProps.OkayCellsize=(meancellradius >= initval.Screen.CellsizeMinRadius);
    GeneralCellProps.OkayCellmorph=(stdcellradius <= initval.Screen.CellsizeMaxStd); 
    if isfield(initval.Screen, 'CellCircularity')
        GeneralCellProps.OkayCellCircularity=(cell_circularity <= initval.Screen.CellCircularity); 
    else
        GeneralCellProps.OkayCellCircularity=1;
    end
    Summary.cellsizeOK=Summary.cellsizeOK+...
        GeneralCellProps.OkayCellsize*GeneralCellProps.OkayCellmorph*GeneralCellProps.OkayCellCircularity;    
    %collect:
    CellSizeVals=[meancellradius stdcellradius GeneralCellProps.OkayCellsize];     
    summary_vals.meancellradius(jj)=meancellradius;
    summary_vals.stdcellradius(jj)=stdcellradius;
    summary_vals.cell_circularity(jj)=cell_circularity;
    
    %% Filter2: Chromosome size and shape
    averageradius_chro = mean(Chromosome.PolarContourMax);
    stdradius_chro = std(Chromosome.PolarContourMax);
    chro_circularity=Chromosome.TotalMaxPeakLength/averageradius_chro/(2*pi);
    %evaluate:
    GeneralCellProps.Okayradiusmean=(averageradius_chro>=initval.Screen.ChromosomeSizeMinRadius);
    GeneralCellProps.Okayradiusstd=(stdradius_chro<=initval.Screen.ChromosomeSizeMaxStd);    
    if isfield(initval.Screen, 'ChroCircularity')
        GeneralCellProps.OkayChroCircularity=(chro_circularity <= initval.Screen.ChroCircularity); 
    else
        GeneralCellProps.OkayChroCircularity=1;
    end
    Summary.chromosomesizeOK=Summary.chromosomesizeOK+...
        GeneralCellProps.Okayradiusmean*GeneralCellProps.Okayradiusstd*GeneralCellProps.OkayChroCircularity;  
    
    %collect:
    ChroSizeVals=[averageradius_chro stdradius_chro GeneralCellProps.Okayradiusmean];        
    summary_vals.averageradius_chro(jj)=averageradius_chro;
    summary_vals.stdradius_chro(jj)=stdradius_chro;
    summary_vals.chro_circularity(jj)=chro_circularity;    
    
    %% Filter3: Is the center at the dark region
    centerintensity = Chromosome.picture(uint8(Chromosome.yCOM),uint8(Chromosome.xCOM));
    maxintensity = max(Chromosome.picture(:));
    threshold = initval.Screen.ChromosomeDonutHoleDepthMin;  %Widefield    
    GeneralCellProps.Okaycenter = (centerintensity<threshold*maxintensity);% Average intensity
    Summary.donutOK=Summary.donutOK+GeneralCellProps.Okaycenter;    
    DonutNess=[centerintensity maxintensity centerintensity/maxintensity GeneralCellProps.Okaycenter];
    summary_vals.donutness(jj)=centerintensity/maxintensity;
    
    %% Filter4: Is the angle ok?
    %MinExpectAngle=GetMinangle(initval); set init in degrees
    MinMeasAngle=mod(abs(c4.spotPosAngle-c3.spotPosAngle)*180/pi,360);
    if MinMeasAngle>180, MinMeasAngle=360-MinMeasAngle; end     
    GeneralCellProps.OkayTerOriAngle=(MinMeasAngle>=initval.Screen.MinimalLabelAngle);
    Summary.teroriangleOK=Summary.teroriangleOK+...
        GeneralCellProps.OkayTerOriAngle;
    summary_vals.ori_ter_angle(jj)=MinMeasAngle;
    
    %% Filter5: Is the distance between Center and locus reasonable?
    Distoricent = sqrt((c3.spotY-Chromosome.yCOM)^2+(c3.spotX-Chromosome.xCOM)^2);
    Disttercent = sqrt((c4.spotY-Chromosome.yCOM)^2+(c4.spotX-Chromosome.xCOM)^2);
    GeneralCellProps.OkayTerCenterDist = ( (Disttercent/averageradius_chro) >= initval.Screen.TerOriMinDistFromCenter);
    GeneralCellProps.OkayOriCenterDist = ((Distoricent/averageradius_chro) >=  initval.Screen.TerOriMinDistFromCenter);
    Summary.terorilociiOK=Summary.terorilociiOK+...
        GeneralCellProps.OkayTerCenterDist*GeneralCellProps.OkayOriCenterDist;
    summary_vals.distance_ori_center(jj)=Distoricent/averageradius_chro*100;
    summary_vals.distance_ter_center(jj)=Disttercent/averageradius_chro*100;
    
    %% Filter6: is the main spot strong enough (i.e., can we reject multi-spot?
    
    GeneralCellProps.OkayOriStrength=(c4_xt.mainpeakfraction>initval.Screen.OriSingleSpotFraction);
    summary_vals.mainpeakfraction(jj)=c4_xt.mainpeakfraction;
    Summary.spotok=Summary.spotok+1.0*GeneralCellProps.OkayOriStrength;
    if 0
        subplot(2,2,1); pcolor(c3_pic); shading flat; axis equal; title('c3');
        subplot(2,2,2); pcolor(c4_pic); shading flat; axis equal; title('c4');
        subplot(2,2,3); pcolor(c3_pic_xt); shading flat; axis equal; title('c3-MAX');
        subplot(2,2,4); pcolor(c4_pic_xt); shading flat; axis equal; title('c4-MAX');
        text(3,3,num2str(1.0*GeneralCellProps.OkayOriStrength), 'Color','w');
        [~]=ginput(1);
    end
    %% Product of all filters
    GeneralCellProps.Okayproduct = GeneralCellProps.OkayCellsize*GeneralCellProps.OkayCellmorph*GeneralCellProps.Okayradiusmean...
        *GeneralCellProps.Okayradiusstd*GeneralCellProps.Okaycenter*GeneralCellProps.OkayTerOriAngle*GeneralCellProps.OkayTerCenterDist...
        *GeneralCellProps.OkayOriCenterDist*GeneralCellProps.OkayOriStrength*GeneralCellProps.OkayByUser;
    
    if initval.PassAllCells==1, GeneralCellProps.Okayproduct=GeneralCellProps.OkayByUser; end
        
    GeneralCellProps=orderfields(GeneralCellProps);
    save(strcat(MatFilePath,CellName),'GeneralCellProps', '-append');
    SummaryTable=[SummaryTable;... 
                  [CellSpecs CellSizeVals ChroSizeVals DonutNess GeneralCellProps.OkayByUser GeneralCellProps.Okayproduct ]];
    
              I=imread(strcat(strcat(initval.resultpath,'CellImages_All',initval.DirSep,''),'Cell',cellno,'.jpg')); 
    if actions.do_the_cell_jpg
        if GeneralCellProps.Okayproduct > 0
            Summary.totalok=Summary.totalok+1;         
            imwrite(I,strcat(strcat(imoutdir,'Cell',cellno,'.jpg')));  
        else
            imwrite(I,strcat(strcat(imoutdir_rej,'Cell',cellno,'.jpg'))); 
        end
    end
    end
end

if ~dryrun
%post-processing 
CellProps_Av=nanmean(SummaryTable);
CellProps_Std=nanstd(SummaryTable);
CellProps_Av(1:2)=[-1 -10000];
CellProps_Std(1:2)=[-2 -10000];

SummaryTable=[CellProps_Av;
            CellProps_Std;
            SummaryTable];

%saving
ColNames=[{'index'} , {'label'}, {'cellradius_mean'} ,{'cellradius_std'}, {'Cellsize OK'},...
           {'chrom_averageradius'}, {'chrom_stdradius'}, {'Okayradius'}, ...
           {'I_COMcenter'} , {'I_max'}, {'Ratio'} , {'DonutOK'}, {'ManualUserOK'},{'CellApproved'}];
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A015_Cell_ScreenReport.xlsx'),ColNames,'Sheet1','A1');
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A015_Cell_ScreenReport.xlsx'),SummaryTable,'Sheet1','A2');
    save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A015_Cell_ScreenReport.mat'),'SummaryTable'); 

    
    
%diagnose plots 
close all;
figure('Units','normalized','Position',[0 0 1 1]);
show_distribution_cutoffs(summary_vals.meancellradius,initval.Screen.CellsizeMinRadius,3,5,1,'cell:R_a_v');
show_distribution_cutoffs(summary_vals.stdcellradius,initval.Screen.CellsizeMaxStd,3,5,2,'cell:R_s_t_d');
show_distribution_cutoffs(summary_vals.cell_circularity,initval.Screen.CellCircularity,3,5,3,'cell:circularity');
show_distribution_cutoffs(summary_vals.averageradius_chro,initval.Screen.ChromosomeSizeMinRadius,3,5,4,'chro:R_a_v');
show_distribution_cutoffs(summary_vals.stdradius_chro,initval.Screen.ChromosomeSizeMaxStd,3,5,5,'chro:R_s_t_d');
show_distribution_cutoffs(summary_vals.chro_circularity,initval.Screen.ChroCircularity,3,5,6,'chro:circularity');
show_distribution_cutoffs(summary_vals.distance_ter_center,100*initval.Screen.TerOriMinDistFromCenter,3,5,7,'labels:radialpos_t_e_r');
show_distribution_cutoffs(summary_vals.distance_ori_center,100*initval.Screen.TerOriMinDistFromCenter,3,5,8,'labels:radialpos_o_r_i');
show_distribution_cutoffs( summary_vals.ori_ter_angle,initval.Screen.MinimalLabelAngle,3,5,9,'labels:ori-ter-angle');
show_distribution_cutoffs(summary_vals.donutness,initval.Screen.ChromosomeDonutHoleDepthMin,3,5,10,'chro:donutdepth');
show_distribution_cutoffs(summary_vals.mainpeakfraction,initval.Screen.OriSingleSpotFraction,3,5,11,'labels:spot fraction');

%saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A015_ScreeningReport_overview.jpg'),'jpg');
% 
%summary plot
subplot(3,2,6);
Summary;
SummaryLegend=([{'1.initialcount'},{'2.cellshapeOK'},{'3.chromosomeshapOK'},{'4.donutdepthOK'},...
               {'5.teroriangleOK'},{'6.terorilociiOK'},{'7.userok'},{'8.SpotOK'},{'9.PRODUCTOK'}]);
SummaryVals=[   Summary.initialcount        Summary.cellsizeOK...
                Summary.chromosomesizeOK   Summary.donutOK...
                Summary.teroriangleOK      Summary.terorilociiOK...
                Summary.userok             Summary.spotok       Summary.totalok];
            
disp('done');
LS=length(SummaryVals);
axz=1:LS;
for ii=1:LS
    bar(axz(ii),SummaryVals(ii)); hold on;    
end
legend(SummaryLegend,'Location','EastOutside');
title(Replace_underscores(strcat('Screening Report of: ',initval.expi)),'FontSize', 8);
disp('(save as jpg and svg....');
saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A015_ScreeningReport.jpg'),'jpg');
pause(0.5);
plot2svg(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A015_ScreeningReport.svg'), gcf);

               
end
disp('done');




function MinExpectAngle=GetMinangle(initval);
if initval.StopLabelpos>initval.StartLabelpos  %ori-start-stop
    a1=(initval.StopLabelpos-initval.StartLabelpos)/100*360; 
    a2=(initval.StartLabelpos+(100-initval.StopLabelpos))/100*360;
else  %start-ori-stop
    a1=((100-initval.StartLabelpos)+initval.StopLabelpos)/100*360; 
    a2=(initval.StartLabelpos-initval.StopLabelpos)/100*360;
    
end
MinExpectAngle=min([a1 a2]);

function UserRejectedCellsList=GetRejectedCellNumbers(initval)
%1 old format, for numeric only
RejectionsFile=strcat(initval.resultpath,'ManuallyRejectedCells');
if exist(strcat(RejectionsFile,'.xlsx'))==2
    [cellnumbers,cellstrings]=xlsread(strcat(RejectionsFile,'.xlsx'));
    if length(cellnumbers)>1 %numeric
        UserRejectedCellsList=cellnumbers;
    else
        LL=length(cellstrings);
        UserRejectedCellsList=zeros(LL,1);
        for kk=1:LL
            cellno=char(cellstrings{kk});
            UserRejectedCellsList(kk)=BuildNumericCellLabel(cellno);
        end
    end
else %check old format
    if exist(strcat(RejectionsFile,'.txt'))==2
        UserRejectedCellsList=dlmread(strcat(RejectionsFile,'.txt'));
    else
        UserRejectedCellsList=[];
    end
end;
dum=1;


function show_distribution_cutoffs(histvals,cutoffs,rws,cols,idx,titl);
nbins=5+ceil(length(histvals)/10);
[this_hist,binax]=hist(histvals,nbins);
subplot(rws,cols,idx);
bar(binax,this_hist); hold on
for ii=1:length(cutoffs)
    markerx=[cutoffs(ii) cutoffs(ii)];
    markery=[0 max(this_hist)];
    plot(markerx,markery, '--', 'LineWidth',2);
end
axis tight;
title(titl, 'FontSize', 6);

