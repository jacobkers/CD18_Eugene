function A013_WF_PerCell_AnalyzeSpotsStandAlone
%JWJK_A:-----------------------------------------------------
%Spot analysis per channel
%Description: different color channels and the widefield channel are
%allocated to the proper analysis functions, depending on the type of
%pattern
%Input: various
%Output: various spot data .mat files
%References: Jacob Kerssemakers, Cees Dekker Lab, Delft
%:JWJK_A------------------------------------------------
close all;
%runtime options
actions.showandsavepics=1;
batchrunindex=102;
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
disp(initval.expi);
%if isdir(initval.resultpath), rmdir(initval.resultpath,'s');  end
%mkdir(initval.resultpath);
mkdir(strcat(initval.resultpath,'CellImages_All',initval.DirSep));
MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
LC=length(initval.Cell_Labels);
AllCellBasicGeometries=[];

for ii=1:LC
    disp(strcat('Analyzing Color channels..,', num2str(LC-ii+1), 'cells to go'));
    %cellno=num2str(initval.cell_labelnos(ii)); 
    cellno=char(initval.Cell_Labels{ii});    
    NumCellLabel=BuildNumericCellLabel(cellno);
    CellSpecs=[ii NumCellLabel];
       
    extraedge=1;
    switch initval.cellmaskname  %SuperUglyHack
        case 'cellmask'; 
        load(strcat(initval.maindatapath,'ma_',cellno,'.mat'), initval.cellmaskname);
        cellmask=bwmorph(cellmask,'dilate',extraedge);
        case 'cellma'
          load(strcat(initval.maindatapath,'ma_',cellno,'.mat'), 'cellma');
          cellmask=bwmorph(cellma,'dilate',extraedge);
    end   
    edge_pic = bwmorph(cellmask,'remove');
%% 
   %------------------------------------------------------------
   %3)The respective ter and ori images ('cfp' and 'rfp') are loaded.  
   % For each, a single spot center and corresponding spot properties 
   % (position, intensity) is tracked following 
   % Llorente-Garcia / Reyes et al.
   % Possible leakage of the chromosome channel in the spot detection
   % channel is background-removed by succesive 'spot peeling'.
   
%% Cfp------------------------------------------------------------------ 
    if 0
        load(strcat(initval.maindatapath,'c1_',cellno,'.mat'),'cellc1');
        Cfp_pic=GetWorkpicFromStack(cellc1,'FocalPlane');
    else
        if strcmp(initval.searchlabel,'c')
            load(strcat(initval.maindatapath,initval.searchlabel,'1_',cellno,'.mat'),'cellc1');
            Cfp_pic=GetWorkpicFromStack(cellc1,'FocalPlane');
         end
         if strcmp(initval.searchlabel,'rc')
            load(strcat(initval.maindatapath,initval.searchlabel,'1_',cellno,'.mat'),'rcellc1');
            Cfp_pic=GetWorkpicFromStack(rcellc1,'FocalPlane');
         end  
    end
    Cfp_pic=Cfp_pic-median(Cfp_pic(:));    
    Cfp_pic=Cfp_pic.*cellmask;
    
    Cfp=Get_MultiSpotProps(Cfp_pic,initval);    
    %% Rfp------------------------------------------------------------------
    if 0
        load(strcat(initval.maindatapath,'c3_',cellno,'.mat'),'cellc3');
        Rfp_pic=GetWorkpicFromStack(cellc3,'FocalPlane');
    else
        if strcmp(initval.searchlabel,'c')
            load(strcat(initval.maindatapath,initval.searchlabel,'3_',cellno,'.mat'),'cellc3');
             Rfp_pic=GetWorkpicFromStack(cellc3,'FocalPlane');
         end
         if strcmp(initval.searchlabel,'rc')
            load(strcat(initval.maindatapath,initval.searchlabel,'3_',cellno,'.mat'),'rcellc3');
             Rfp_pic=GetWorkpicFromStack(rcellc3,'FocalPlane');
         end   
    end
    Rfp_pic=Rfp_pic-median(Rfp_pic(:));
    Rfp_pic=Rfp_pic.*cellmask;
    
    Rfp=Get_MultiSpotProps(Rfp_pic,initval);
    
   if actions.showandsavepics 
           pcolor(1.0*edge_pic); shading flat, colormap bone; hold on; 
           %plot(pairs1,pairs2,'w-'); hold on;    
           plot(Cfp.spotX,Cfp.spotY, 'bo','MarkerSize',10, 'MarkerFaceColor','b'); hold on;
           plot(Rfp.spotX,Rfp.spotY, 'ro','MarkerSize',10,'MarkerFaceColor','r'); hold on;

           title(strcat('Cell', num2str(cellno,'% 3.0f')));
           pause(0.01);

            saveas(gcf,strcat(initval.resultpath,'CellImages_All',initval.DirSep,'Cell', num2str(cellno,'% 3.0f'),'.jpg')); 
            pause(0.3); 
            %[~]=ginput(1); 
            close(gcf);
    end
    initval.ResultName=strcat(MatFilePath,strcat('ResultsOfCell',cellno,'_Spots.mat'));
    save(initval.ResultName, 'Cfp','Rfp','Cfp_pic','Rfp_pic');
end

function pic=GetWorkpicFromStack(stack,WorkPicOption);
    stack=double(stack);
    switch WorkPicOption
        case 'MeanProject',pic=squeeze(mean(stack,3));
        case 'FocalPlane', 
        %get main plane  (assuming largely planar features)
        stcurve=squeeze(std(std(stack)));
        [~,MainPlane]=max(stcurve);
        pic=double(stack(:,:,MainPlane));  
        case 'First', pic=double(stack(:,:,1));
    end
       
 function spot=Get_MultiSpotProps(pic,initval);
    Psf=initval.Psf_est;
    %get a list of tentative spot positions
    AllSpotProps=PeelblobsFromImage(pic,Psf,0.8,0.1,0);
    [Nsp,~]=size(AllSpotProps);
    for ii=1:Nsp %refine spots
    
    CellSimProps.roilox=1;          %lower x coordinate of Roi
    CellSimProps.roiloy=1;           %lower y coordinate of Roi
    CellSimProps.absx=AllSpotProps(ii,3);    %first estimate
    CellSimProps.absy=AllSpotProps(ii,4);    %first estimate     
    [spot_fit,spot_est,spotim, bckim]=SpotsBy_1x2DGaussFixWidth_BackgroundBy_localROI_iterative(pic,CellSimProps,Psf,0);             
    spot.spotY(ii)=spot_fit.y0;
    spot.spotX(ii)=spot_fit.x0; 
    spot.spotContent(ii)=spot_fit.N0;  
    end
    
        
function PP= TwoDGaussNormPeak(im,x0,y0,psf)
%This is the equation for a 2D gaussian
[r,c]=size(im);
[XX,YY]=meshgrid(1:c,1:r);
RR=((XX-x0).^2+(YY-y0).^2).^0.5;  %distance of allpixels to clickpoint
PP =exp (-(RR).^2./(2*psf.^2));

       
