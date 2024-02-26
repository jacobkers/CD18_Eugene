function A010_WF_PerCell_AnalyzeAllColors(batchrunindex)
%JWJK_A:-------------------------------------------------------------------
%Title:Position and content analysis of fluorescence colour channels

%General summary : The various color channels
%(chromatin, c3 and c4 labeled positions) are analysed for positions and
%fluorescence counts. The circular (angular) coordinate system to describe
%the chromatin donut shape is set up.
%Approach
%Chromosome channel: Chromosome channel image is loaded. This can be the focal 
%plane of a stack or its maximum projection. This image is masked with the
%brightfield-based mask to remove chromosome intensity from adjacent cells. 
%Next, we perform 'QI' style center-tracking on the ring structure to get 
%the center and a radial map of the chromosome pattern similar to ref [1]. 
%Each radial section yields a radial profile, which is subjected to peak 
%analysis (position,content, width etc. 
%
%2)%c3/c4 channels are presumed to contain one spot each, these are
%located following Llorente-Garcia / Reyes et al [2]. We use a squared
%region, and for the 'putative spot' area pick a random
 %coordinate within some pixels distance of the spot.
%Input: .mat data from precursor analysis, including a mask encompassing 
%the cell contour is used to exclude signal from adjacent cells.
%
%Output: 
%1)data is saved to a structure 'Chromosome' in .mat for further use; 
%images are saved showing the result of the channel analysis.  
%2) an Excel table listing basic geometries per cell
%
%References
%[1] M.T.J. van Loenhout, J. Kerssemakers , I. De Vlaminck, C. Dekker
%Non-bias-limited tracking of spherical particles, enabling nanometer 
%resolution at low magnification, Biophys. J. 102, Issue 10, 2362 (2012)
%[2] Llorente-Garcia I. et al,  Biochim. Biophys. Acta. 2014;1837:811–824
%:JWJK_A-------------------------------------------------------------------

close all;
%runtime options
actions.workfullstackdiagnosis=0;  %works only if also done in cropcode
actions.showandsavepics=1;
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
disp(initval.expi);
%if isdir(initval.resultpath), rmdir(initval.resultpath,'s');  end
%mkdir(initval.resultpath);
mkdir(strcat(initval.resultpath,'CellImages_All',initval.DirSep));
MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
mkdir(MatFilePath);

LC=length(initval.Cell_Labels);
AllCellBasicGeometries=[];



for ii=1:LC
    disp(strcat('Analyzing Color channels:cell',num2str(ii),'and', num2str(LC-ii+1), 'cells to go'));
    %cellno=num2str(initval.cell_labelnos(ii)); 
    cellno=char(initval.Cell_Labels{ii});    
    NumCellLabel=BuildNumericCellLabel(cellno);
    CellSpecs=[ii NumCellLabel];
    
    %load and identify color channels:
    [channel_stack,channel_stack_extra]=get_channel_stack_from_cropcode(cellno,initval,actions); 
    [cell_pic, cell_mask, c2_pic, c3_pic,c4_pic]=get_channel_ID(channel_stack, initval);
    [~, ~,                c2_pic_xt, c3_pic_xt,c4_pic_xt]=get_channel_ID(channel_stack_extra, initval); 
    
    %some further refinement:
    celledge_pic = bwmorph(cell_mask,'remove');    
    c2_pic=(c2_pic-min(c2_pic(:))).*cell_mask;       
    QI_chro=TrackXY_by_QI_Init_Chro(c2_pic,initval);  
    Chromosome=ChromosomeRadialSampling_WideField(c2_pic.*cell_mask,QI_chro, 0); 
    Chromosome.picture=c2_pic;
   %------------------------------------------------------------
   %3)The respective ter and ori images ('c4' and 'c3') are loaded.  
   % For each, a single spot center and corresponding spot properties 
   % (position, intensity) is tracked following 
   % Llorente-Garcia / Reyes et al.
   % Possible leakage of the chromosome channel in the spot detection
   % channel is background-removed by succesive 'spot peeling'.
   
    % c3 spots
    c3_pic=c3_pic-median(c3_pic(:));
    c3_pic=c3_pic.*cell_mask;    
    c3=Get_SpotProps(c3_pic,Chromosome,initval);
    % c4 spots:
    c4_pic=c4_pic-median(c4_pic(:));    
    c4_pic=c4_pic.*cell_mask;    
    c4=Get_SpotProps(c4_pic,Chromosome,initval);    
    % c3_xt spots
    c3_pic_xt=c3_pic_xt-median(c3_pic_xt(:));
    c3_pic_xt=c3_pic_xt.*cell_mask;    
    c3_xt=Get_SpotProps(c3_pic_xt,Chromosome,initval);
    %c4_xt spots:
    c4_pic_xt=c4_pic_xt-median(c4_pic_xt(:));    
    c4_pic_xt=c4_pic_xt.*cell_mask;    
    c4_xt=Get_SpotProps(c4_pic_xt,Chromosome,initval);    
    
    
    
    %% 
    
    % 'Use this section for a Quicksheet'
   %------------------------------------------------------------
   %4)Next, for each sampled chromosome radial section, the corresponding 
   % radial position of the cell wall is processed from the brightfiled
   % edge image. After that, som e basic cell properties are collected
% 'End of Quicksheet section'
     Chromosome=LinkCellWallToChromosomeEdge(Chromosome,celledge_pic);
     
     %diagnosis (hidden optional)
     if actions.workfullstackdiagnosis
         if ii==1, Z_info=struct('Plane_av_raw',[]);end
         Z_info=JK_diagnosis_GetVerticalInfo(Chromosome,chro_stack_raw,chro_stack_decon,c2_pic,initval,cellno,Z_info,ii,LC);
     end

     ThisCellBasicGeometry=Get_BasicCellGeometry(Chromosome);
     
     %Get some general pattern properties
     [fluo,modelpic]=Processing_Fluorescence_SimplePatternAnalysis(Chromosome.picture);
     
     AllCellBasicGeometries=[AllCellBasicGeometries; ...
                             CellSpecs ThisCellBasicGeometry ...
                             fluo.content_total c3.all_image_content c4.all_image_content];

     
     %Get some general pattern properties
     [fluo,modelpic]=Processing_Fluorescence_SimplePatternAnalysis(Chromosome.picture);
                         
    dum=1;
    if 0
    figure(1);
    subplot(3,3,1); pcolor(double(celledge_pic)); colormap hot; shading flat;
    subplot(3,3,2); 
        pcolor(double(chro_pic)); colormap hot; shading flat; hold on;
        plot(Chromosome.xCOM, Chromosome.yCOM, 'yx','MarkerSize', 20); hold off;
    subplot(3,3,3);  
        pcolor(double(c4_pic)); colormap hot; shading flat; hold on;
        plot(c4.spotX, c4.spotY, 'bx','MarkerSize', 20); hold off;
        
    subplot(3,3,4); 
        pcolor(double(c3_pic)); colormap hot; shading flat; hold on;
        plot(c3.spotX, c3.spotY, 'rx','MarkerSize', 20); hold off; 
    
    subplot(3,3,5); pcolor(chro_pic); colormap hot; shading flat; hold on;
    
    initval.ResultName=strcat(MatFilePath,strcat('ResultsOfCell',cellno,'.mat'));
    save(initval.ResultName, 'celledge_pic','cell_pic','cellmask', 'chro_pic','c4_pic','c3_pic');
    [~]=ginput(1);
    end
 
    
   if actions.showandsavepics 
            set(figure(2), 'visible','off')            
            subplot(2,4,1);
                pcolor(double(cell_mask)); colormap hot; shading flat; hold on;
                plot(Chromosome.xCOM, Chromosome.yCOM, 'yx','MarkerSize', 20); hold off;
                title('mask');
                axis equal; axis tight, axis off;
            subplot(2,4,2);
                pcolor(double(c2_pic)); colormap hot; shading flat; hold on;
                plot(Chromosome.xCOM, Chromosome.yCOM, 'yx','MarkerSize', 20); hold off;
                title('c2/DNA');
                axis equal; axis tight, axis off;
            subplot(2,4,5);
                pcolor(double(c3_pic)); colormap hot; shading flat; hold on;
                plot(c3.spotX, c3.spotY, 'rx','MarkerSize', 15); hold off; 
                axis equal; axis tight, axis off;
                title('c3')
            subplot(2,4,6); 
                pcolor(double(c4_pic)); colormap hot; shading flat; hold on;
                plot(c4.spotX, c4.spotY, 'bx','MarkerSize', 15); hold off;
                axis equal; axis tight, axis off; 
                title('c4')
           subplot(1,2,2);
           NoGaps=Chromosome.ValidIdx;
           CX=Chromosome.CartesianContourEdge_X;
           CY=Chromosome.CartesianContourEdge_Y;
           MX=Chromosome.CartesianContourMax_X;
           MY=Chromosome.CartesianContourMax_Y;
           BX=Chromosome.CartesianCellwallX;
           BY=Chromosome.CartesianCellwallY;
           pcolor(Chromosome.picture); colormap bone; shading flat; hold on;
           plot(Chromosome.xCOM,Chromosome.yCOM,'yx-', 'MarkerSize',10); hold on;
           plot(CX(NoGaps),CY(NoGaps),'y-','LineWidth',1); hold on;
           plot(MX(NoGaps),MY(NoGaps),'y-','LineWidth',1); hold on;
           plot(BX,BY,'w-', 'LineWidth',1); hold on;
           NBX=BX(Chromosome.CellEdgeNearestIndex);
           NBY=BY(Chromosome.CellEdgeNearestIndex);
           pairs1=[CX(NoGaps)';NBX(NoGaps)']; 
           pairs2=[CY(NoGaps)';NBY(NoGaps)'];
           %plot(pairs1,pairs2,'w-'); hold on;    
           plot(c4.spotX,c4.spotY, 'bo','MarkerSize',6, 'MarkerFaceColor','b'); hold on;
           plot(c3.spotX,c3.spotY, 'ro','MarkerSize',6,'MarkerFaceColor','r'); hold on;

           title(Replace_underscores(strcat('Cell', num2str(cellno,'% 3.0f'))));
           axis equal; axis tight; 
           pause(0.01);

            saveas(gcf,strcat(initval.resultpath,'CellImages_All',initval.DirSep,'Cell', num2str(cellno,'% 3.0f'),'.jpg')); 
            pause(0.01); 
            %[~]=ginput(1); 
            close(gcf);
    end
    initval.ResultName=strcat(MatFilePath,strcat('ResultsOfCell',cellno,'.mat'));
    save(initval.ResultName, 'Chromosome','c4','c3','c4_xt','c3_xt',...
                            'celledge_pic','cell_pic','cell_mask', ...
                            'c2_pic','c4_pic','c3_pic','c4_pic_xt','c3_pic_xt');
end



%post-processing 

CellProps_Av=nanmean(AllCellBasicGeometries);
CellProps_Std=nanstd(AllCellBasicGeometries);
CellProps_Av(1:2)=[-1 -10000];
CellProps_Std(1:2)=[-2 -10000];

AllCellBasicGeometries=[CellProps_Av;
                        CellProps_Std;
                        AllCellBasicGeometries];
         
ColNames=[{'index'} , {'label'}, ...
        {'chromosome_peak_length'}, ...
        {'chromosome_outercontour_length'}, ...
        {'chrom_peak_averageradius'},{'chrom_peak_stdradius'},...
        {'chrom_peak_maxradius'},{'chrom_peak_minradius'},...
        {'chrom_outercontour_averageradius'},{'chrom_outercontour_stdradius'},...
        {'chrom_outercontour_maxradius'},{'chrom_outercontour_minradius'},...
        {'chrom_radiusofgyration'},...
        {'chrom_averageFWHM'},{'chrom_stdFWHM'},...
        {'chrom_maxFWHM'},{'chrom_minFWHM'},...
        {'cell wall length'},{'cellwall_averageradius'},{'cellwall_stdradius'},...
        {'cellwall_maxradius'},{'cellwall_minradius'}, ...
        {'rawcounts_chro'}, {'rawcounts_c3'},{'rawcounts_c4'}];
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A010_Cell_BasicGeometries.xlsx'),ColNames,'Sheet1','A1');
pause(1);
xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A010_Cell_BasicGeometries.xlsx'),AllCellBasicGeometries ,'Sheet1','A2');
pause(1);
save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A010_Cell_BasicGeometries.mat'),'AllCellBasicGeometries');             



    

     function Chromosome=LinkCellWallToChromosomeEdge(Chromosome,celledgepic);              
       %from chromosome
       NoGaps=find((1.0-Chromosome.Gaps)==1);
       CHangle=Chromosome.AnnularAxis;
       CHx0=Chromosome.xCOM;
       CHy0=Chromosome.yCOM;
       
       
       %from cell edge
       [rr,cc]=size(celledgepic);
       [XX,YY]=meshgrid(1:cc,1:rr);
       sel=find(celledgepic>0);
       BX=XX(sel);
       BY=YY(sel);
       
       CX=Chromosome.CartesianContourEdge_X;
       CY=Chromosome.CartesianContourEdge_Y;
       
       
       %% 1 get radial cell edge distance to chromosome center, sorted 0-2*pi
       BangleC1=atan2((BY-CHy0),(BX-CHx0));  %angles
       BangleC2=mod(BangleC1,2*pi);
       [BangleC,ix]=sort(BangleC2);            %sorted 0-2pi   
       BradC=((BX-CHx0).^2+(BY-CHy0).^2).^0.5; %radial distances
       BradC=BradC(ix);                        %sorted 0-2pi 
       CellEdgedistFromChomosomeCenter=interp1(BangleC, BradC,CHangle);  %radii mapped on chromosome angles
       Chromosome.RadialCellEdge=CellEdgedistFromChomosomeCenter;
       Chromosome.CartesianCellwallX=Chromosome.xCOM+CellEdgedistFromChomosomeCenter.*cos(CHangle);
       Chromosome.CartesianCellwallY=Chromosome.yCOM+CellEdgedistFromChomosomeCenter.*sin(CHangle);       
       
       %% 2 get the nearest cell wall point per chromosome position
       LC=length(CX);
       CellEdgedistNearest=zeros(LC,1);
       CellEdgedistNearestIndex=zeros(LC,1);
       for cc=1:LC  %for  all cell positions
           cx=CX(cc);
           cy=CY(cc);
           [mindist,ix]=min(((Chromosome.CartesianCellwallX-cx).^2+(Chromosome.CartesianCellwallY-cy).^2).^0.5); %minimum radial distance
           CellEdgeNearestDist(cc)=mindist;   %minimum distance
           CellEdgeNearestIndex(cc)=ix;       %minimum point (index)
            %minimum edge
       end
       Chromosome.CellEdgeNearestDist=CellEdgeNearestDist;
       Chromosome.CellEdgeNearestIndex=CellEdgeNearestIndex;      
       %% 
       

 
    
     
