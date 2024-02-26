function A050_WF_2DClusterAnalysis(batchrunindex)
%JWJK_A:-------------------------------------------------------------------
%Two-dimensional (image-based) cluster analysis
%
%Summary: This function performs 2D cluster analysis on chromatin patterns of 
%cells and stores the result per cell, per cluster.
%
%Approach: we deconstruct the image by means of single-spot, psf-sized
%Gaussians. This is done by subtracting such Gaussians from the image by
%fist picking the brightest image point and subtracting one Gaussian of
%equal peak intensity there, then repeating this action on the residual
%image and continuing to do so until all intensity is covered by the subtracted 
%Gaussians.
%Next, the Gaussians are grouped in clusters. Each cluster consists of a
%small group of Gaussians(typically 1-5) or 'components' that are within one psf 
%distance of one of the others, such that the cluster forms an optically connected shape. 
%Likewise, all components of one cluster are at least one psf away form
%those of another.
%These clusters are then parametrized by their center-of mass position,
%content, radius of gyration and so on. Finally, their angular position with respect to  
% the chromatin geometrical center (as was determined via the 1D density curve analysis)
% is determined and stored.
%
%Input: data in .mat files stored in former analysis steps.
%
%Output: Data is presented as scatter plots and histograms. Saved are 
%Tabular excel files, .mat and summary plots.
%
%:JWJK_A-------------------------------------------------------------------

close all;

justplotit=0;
if nargin<1,batchrunindex=1; end
initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
if initval.UseMeasuredPSFforClusterAnalysis
    Psf_meas=B050_EstimatePSF(batchrunindex)
else
    Psf_meas=initval.Psf_est
end

allframes=length(initval.Cell_Labels);
imoutdir=strcat(initval.resultpath,'CellImages_Clusteranalysis',initval.DirSep);   
%if isdir(imoutdir), rmdir(imoutdir,'s');  end
mkdir(imoutdir);
MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);

all_channel_suffixes=[{'_c3'},{'_c4'},{'_c2'}];
for ci=1:3
goodcount=0;
ClusterMaxNo=0;
channel_suffix=char(all_channel_suffixes{ci});
%analysis loop
AllCellsAllClusters=[];
AllClusters=[];
ClusterGenprops=struct('cell_label',[]);
if ~justplotit
for jj=1:allframes    
    cellno=char(initval.Cell_Labels{jj});     
    CellName=strcat('ResultsOfCell',cellno,'.mat'); 
        load(strcat(MatFilePath,CellName));   
        CellOK=(GeneralCellProps.Okayproduct); 
        if (CellOK|initval.PassAllCells) 
        disp(strcat(initval.expi,'_',channel_suffix,'_',CellName,'ClusterAnalysis..', num2str(allframes-jj+1), 'cells to go'));
        goodcount=goodcount+1;        
        cellno=char(initval.Cell_Labels{jj});    
        NumCellLabel=BuildNumericCellLabel(cellno);
        CellSpecs=[goodcount NumCellLabel]; 
        %work the clusters:
        switch channel_suffix
            case '_c2', cluster_pic=c2_pic;  %chromosome
            case '_c3', cluster_pic=c3_pic;  %c3
            case '_c4', cluster_pic=c4_pic;  %c4
        end
        cluster_pic=Remove_Background(cluster_pic,'Min');        
        AllSpotProps=PeelblobsFromImage(cluster_pic,Psf_meas,0.98,0.03, 0); 
        AllSpotProps=F006_CleanSpots(AllSpotProps,cluster_pic,Psf_meas);      
        [rr,cc]=size(AllSpotProps);        
            if rr>1
            x=AllSpotProps(:,3);
            y=AllSpotProps(:,4);
            ClusterBasics=Find_Clusters(AllSpotProps, initval);
            [~,CL]=size(ClusterBasics);
           
            %Build single cluster contours
            [Clusters,ClusterPropsRow,ThisCellClusterTable,AllContours]=GetClusterDetails(ClusterBasics,Chromosome,Aligned,cluster_pic,cellno);
           
            %Get some general properties per cell 
            clus_pic=AllContours.ReconstructIm;
            ClusterGenprops=Get_ClusterGenprops(ClusterGenprops,ThisCellClusterTable,cluster_pic,clus_pic,CellSpecs);

            [~,clusterno]=size(Clusters);
            ClusterMaxNo=max([ClusterMaxNo clusterno]); %needed for table header

            %add to full per-cell-table, elongate if necessary
            [rr,cc]=size(AllCellsAllClusters);
            newrowentry=[CellSpecs clusterno ClusterPropsRow]; le=length(newrowentry);
            cc_new=max([cc,le]);        
            newtable=zeros(rr+1,cc_new); 
            newtable(1:rr,1:cc)=AllCellsAllClusters;
            newtable(rr+1,1:le)=newrowentry;
            AllCellsAllClusters=newtable;

            %add to cluster-sorted table:
            AllClusters=[AllClusters; ...
                        [repmat(CellSpecs,clusterno,1)...
                         ThisCellClusterTable]];   
            %save cluster info for this cell: 
            switch channel_suffix
                case '_c2' 
                    Clusters_c2=Clusters;
                    AllContours_c2=AllContours;
                    save(strcat(MatFilePath,CellName),'Clusters_c2','AllContours_c2','-append');  %chromosome
                case '_c3' 
                    Clusters_c3=Clusters;
                    AllContours_c3=AllContours;
                    save(strcat(MatFilePath,CellName),'Clusters_c3','AllContours_c3','-append');  %chromosome%c3
                case '_c4', 
                    Clusters_c4=Clusters;
                    AllContours_c4=AllContours;
                    save(strcat(MatFilePath,CellName),'Clusters_c4','AllContours_c4','-append');  %chromosome%c4
            end
             
            end  
        end
end

%% wrapping up for saving: cluster report
    [rr,cc]=size(AllCellsAllClusters);
    CellProps_Av=zeros(1,cc);
    CellProps_Std=zeros(1,cc);
    for ii=1:cc
        bufcol=AllCellsAllClusters(:,ii); bufcol=bufcol(bufcol>0);
        CellProps_Av(ii)=nanmean(bufcol);
        CellProps_Std(ii)=nanstd(bufcol);
    end
    CellProps_Av(1:2)=[-1 -10000];
    CellProps_Std(1:2)=[-2 -10000];

    AllCellsAllClusters=[CellProps_Av;
                    CellProps_Std;
                    AllCellsAllClusters];

    [Header,Headershort]=Build_cluster_Excel_Header(ClusterMaxNo);
        
    [rr,cc]=size(AllCellsAllClusters);
    erasersheet=NaN*zeros(2E3,2*cc);
    xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_Cell_ClusterReport',channel_suffix,'.xlsx'),Header,'Sheet1','A1');
    xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_Cell_ClusterReport',channel_suffix,'.xlsx'),erasersheet,'Sheet1','A3');
    xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_Cell_ClusterReport',channel_suffix,'.xlsx'),AllCellsAllClusters,'Sheet1','A3');
    save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_Cell_ClusterReport',channel_suffix,'.mat'),'AllCellsAllClusters');

    [rr2,cc2]=size(AllClusters);
    erasersheet2=NaN*zeros(10E4,2*cc2);
    xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_ClusterColumnReport',channel_suffix,'.xlsx'),Headershort,'Sheet1','A1');
    xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_ClusterColumnReport',channel_suffix,'.xlsx'),erasersheet2,'Sheet1','A2');
    xlswrite(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_ClusterColumnReport',channel_suffix,'.xlsx'),AllClusters,'Sheet1','A2');
    save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_ClusterColumnReport',channel_suffix,'.mat'),'AllClusters');

    LL=length(ClusterGenprops.clusterno);
    ClusterGenprops.clusterno_av_std=[nanmean(ClusterGenprops.clusterno)...
                                      nanstd(ClusterGenprops.clusterno)...
                                      2*nanstd(ClusterGenprops.clusterno)/LL.^0.5];
    ClusterGenprops.offpercsq_md_std=[nanmedian(ClusterGenprops.off_perc_sq)...
                                      nanstd(ClusterGenprops.off_perc_sq)...
                                      2*nanstd(ClusterGenprops.off_perc_sq)/LL.^0.5];
                                        
    save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_ClusterGeneralprops',channel_suffix,'.mat'),'ClusterGenprops');
    
close all;
end

%% PLOT MENU---------------------------------------------------------------   
%aplot loop
set(figure(1), 'visible','off');
for jj=1:allframes    
    cellno=char(initval.Cell_Labels{jj});     
    CellName=strcat('ResultsOfCell',cellno,'.mat'); 
        load(strcat(MatFilePath,CellName));
        
        CellOK=(GeneralCellProps.Okayproduct);
        %DothisCell=(CellOK|initval.PassAllCells));  
        if (CellOK|initval.PassAllCells) 
            %disp('good cell')
            disp(strcat('Plotting..',CellName,'_',num2str(allframes-jj+1),'cells to go'));
            %set(figure(1), 'visible','on')
            %work the clusters for this channel:

             switch channel_suffix
                    case '_c2'
                        Clusters=Clusters_c2;
                        AllContours=AllContours_c2;
                        cluster_pic=c2_pic;  %chromosome
                        cluster_color='y';
                    case '_c3'
                        Clusters=Clusters_c3;
                        AllContours=AllContours_c3;
                        cluster_pic=c3_pic;  %c3
                        cluster_color='m';
                    case '_c4' 
                        Clusters=Clusters_c4;
                        AllContours=AllContours_c4;
                        cluster_pic=c4_pic;  %c4
                        cluster_color='g';
             end
            plot_pic=cluster_pic-min(cluster_pic(:));
            pcolor(plot_pic'); shading flat, colormap bone; hold on;  
            axis equal; axis tight; axis xy; hold on; 
            PlotProps=Clusters.spotprops;   
            perc=num2str(round(100*PlotProps(end,7)));
            title(strcat('Overlay covering',perc,'percent'));
            
            NoGaps=Chromosome.ValidIdx;      
            %add ridge loop
            MX=Chromosome.CartesianContourMax_X;
            MY=Chromosome.CartesianContourMax_Y;   
            plot(MY(NoGaps),MX(NoGaps),'r-','LineWidth',2); hold on;
            contour(cluster_pic',6,'c','Linewidth',1);
            plot(AllContours.Y,AllContours.X,[cluster_color,'-'],'LineWidth',2); hold off;
            title(strcat('Cell', num2str(cellno,'% 3.0f')));
            pause(0.01);
            saveas(gcf,strcat(imoutdir,'CellClusters', num2str(cellno,'% 3.0f'),channel_suffix,'.jpg')); 
            pause(0.01); 
        end    
end
end

function Clusters=Find_Clusters(SpotProps, initval);
% 'Use this section for a Quicksheet'
%------------------------------------------------------------
    % sort spots by brightness; pick brightest one
    % find near ones for brightmost one from leftover list
    % for the ones found, keep finding new near ones until nochange; 
%------------------------------------------------------------[JK16]
% 'End of Quicksheet section'

 left_ones=SpotProps;
 Clusters=struct([]);
 N_clust=0;
 ClusterTable=[];
 while ~isempty(left_ones);    
    %build a new cluster
    N_clust=N_clust+1;
    Bri=left_ones(:,6); [Brisort,idx]=sort(Bri,'descend');
    thiscluster=left_ones(1,:);
    left_ones=left_ones(2:end,:);  %pick from stock
    [LC,~]=size(thiscluster);
    cluster_growth=1;
    while cluster_growth
        new_ones=[];
        for ii=1:LC  
            %for all existing elements of this cluster, 
            %,find near ones in the leftover stock
            x0=thiscluster(ii,3);
            y0=thiscluster(ii,4);
            leftx=left_ones(:,3);
            lefty=left_ones(:,4);
            rr=((leftx-x0).^2+(lefty-y0).^2).^0.5;
            SeparateWidth=2.5;
            near_ones_idx=find(rr<SeparateWidth*initval.Psf_est); %overlapping
            left_ones_idx=find(rr>=SeparateWidth*initval.Psf_est); 
          
            if ~isempty(left_ones_idx)
                new_ones=[new_ones; left_ones(near_ones_idx,:)]; %add 
                left_ones=left_ones(left_ones_idx,:);  %shrink remaining                 
            end
        end
        if ~isempty(new_ones) 
            thiscluster=[thiscluster; new_ones]; %add new ones to cluster 
        else
            cluster_growth=0;  %...or stop this cluster
        end 
    end
    %add values to cluster struct; go to next
    Clusters(N_clust).spotprops=thiscluster;    
 end


function [Clusters,ThisCellClusterRow,ThisCellClusterTable,AllContours,ReconstructIm]=GetClusterDetails(Clusters,Chromosome, Aligned,orim,cellno)
% 'Use this section for a Quicksheet'
%------------------------------------------------------------
    % this function works trhough all clusters in a cell:
    % builds a single-cluster image and obtains properties for this.
    %input: AllSpotProps, collection of contributing Gaussian spots
    % It also builds a clustercontour(level) using the FWHM lelvel 
    % points (or other levels)
    %ThisCellClusterTable contains sorted  'CL' rows of:
    %order content% area radius density xpos ypos 1Ddistpos 1DBPpos
%------------------------------------------------------------[JK16]
% 'End of Quicksheet section'
 
    [rr,cc]=size(orim);
    [XX,YY]=meshgrid(1:cc,1:rr);    
    [~,CL]=size(Clusters);
    AllContoursX=zeros(50,CL);
    AllContoursY=zeros(50,CL);
    ThisCellClusterTable=zeros(CL,11);
    %contains 'CL' rows of:
    %order content% area radius density xpos ypos 1Ddistpos 1DBPpos
    ReconstructIm=0*orim;   
     for ii=1:CL  %FOR ALL CLUSTERS             
        ThisClusterSpots=Clusters(ii).spotprops;  %components of cluster
        thisclusterim=GetSingleClusterImage(ThisClusterSpots,orim); 
        ReconstructIm=ReconstructIm+thisclusterim;
        C_perc=sum(ThisClusterSpots(:,6));
        psf_used=ThisClusterSpots(1,5);        
        [x_com,y_com,~,~,rad_gyr]=JKD2_IM_calculate2Dmoment_extended(thisclusterim);
        if 0 %strcmp(cellno,'200020')
            pcolor(thisclusterim'); shading flat, colormap bone;
            title(['Cluster' num2str(ii)]);
             axis equal; axis tight; axis xy; hold on;
            [~]=ginput(1);
        end
        [x_1D_BP,x_1D_DS]=Get1DAxisPosition(Chromosome,Aligned,x_com,y_com);                                
        [contourX,contourY,ClusterShapeProps]=GetClusterShapeProps(thisclusterim,orim,x_com,y_com,0.2); 
        
        AllContoursX(:,ii)=contourX;
        AllContoursY(:,ii)=contourY;
        
        
        %add properties to clustertable
        thisclusterprops=[C_perc.... 
                          2*rad_gyr ...                         
                          ClusterShapeProps.EquivDiameter...
                          2*ClusterShapeProps.RgBW... 
                          ClusterShapeProps.Area... 
                          ClusterShapeProps.Density ...
                          x_com... 
                          y_com... 
                          x_1D_DS... 
                          x_1D_BP... 
                          psf_used...
                          ];
        %order content% area radius density xpos ypos 1Ddistpos 1DBPpos
        ThisCellClusterTable(ii,:)=thisclusterprops;
      
        %add properties to 'Cluster' structure
        Clusters(ii).COM_X=x_com;  %[x,y];
        Clusters(ii).COM_Y=y_com;  %[x,y];
        Clusters(ii).COM_BasePerc=x_1D_BP;  %related to ori align
        Clusters(ii).COM_DistPerc=x_1D_DS;  %related to ori align 
        Clusters(ii).C_perc=C_perc;                       
        Clusters(ii).EquivDiameter=ClusterShapeProps.EquivDiameter;
        Clusters(ii).Area=ClusterShapeProps.Area; 
        Clusters(ii).Density=ClusterShapeProps.Density;
        Clusters(ii).psf_used=psf_used; 
     end
     AllContours.X=AllContoursX;
     AllContours.Y=AllContoursY;
     AllContours.ReconstructIm=ReconstructIm;
     
     %Sort cluster table by brightest cluster first, transform  to one row
     %content% radius area density xpos ypos 1Ddistpos 1DBPpos
     [rr,cc]=size(ThisCellClusterTable);
     [~,idx]=sort(ThisCellClusterTable(:,1),'descend');
     ThisCellClusterTable=ThisCellClusterTable(idx,:);
     ThisCellClusterRow=reshape(ThisCellClusterTable',1,rr*cc); %repeat=8;
      
     
     
  function [contourX,contourY,ClusterShapeProps]=GetClusterShapeProps(thisclusterim,orim,xm,ym,cutoff);
      %This function obtains contour points of an spheroid contour, equally
      %spaced, in order of angular revolution around the COM and with a
      %fixed number of points
                BW=0*thisclusterim;
                [rr,cc]=size(thisclusterim);
                [XX,YY]=meshgrid(1:cc,1:rr);
                sel=find(thisclusterim>cutoff*max(orim(:)));
                if length(sel)>1
                BW(sel)=1;   
                ClusterShapeProps = regionprops(BW, 'Area','EquivDiameter');
                ClusterShapeProps.Density=sum(thisclusterim(sel))/length(sel);
                
                [xBW,yBW,~,~,RgBW]=JKD2_IM_calculate2Dmoment_extended(1.0*BW);
                ClusterShapeProps.xBW=xBW; 
                ClusterShapeProps.yBW=yBW; 
                ClusterShapeProps.RgBW=RgBW; 
                
                BWedge=bwmorph(BW,'remove');
                contourX=XX(BWedge);
                contourY=YY(BWedge);
                
                angle=atan2(contourY-ym,contourX-xm);
                [~,idx]=sort(angle);
                contourX=contourX(idx);
                contourY=contourY(idx);
                [contourX,contourY]=B002_EqualizeAlongContour(contourX,contourY,49); 
                %close contour lines by repeating start point at end
                contourX=[contourX ;contourX(1)];
                contourY=[contourY ;contourY(1)];
                
                else
                    contourX=NaN*ones(50,1);
                    contourY=NaN*ones(50,1);      
                    ClusterShapeProps.Area=NaN;
                    ClusterShapeProps.EquivDiameter=NaN;
                    ClusterShapeProps.Density=NaN;
                    ClusterShapeProps.xBW=NaN; 
                    ClusterShapeProps.yBW=NaN; 
                    ClusterShapeProps.RgBW=NaN; 
                end
                
   function thisclusterim=GetSingleClusterImage(ThisClusterSpots,orim);
        %1) build a one-cluster-image
        thisclusterim=0*orim;
         [Nspots,~]=size(ThisClusterSpots);
        for jj=1:Nspots
            Xpos=ThisClusterSpots(jj,3);
            Ypos=ThisClusterSpots(jj,4);
            Psf=ThisClusterSpots(jj,5);
            Peak=ThisClusterSpots(jj,2);
            thisclusterim=thisclusterim+Peak*TwoDGaussNormPeak(orim,Xpos,Ypos,Psf);   
        end
        
  function [GenomicPosition,SpatialPosition]=Get1DAxisPosition(Chromosome,Aligned,xm,ym);
        %find corresponding BP_from-ori and distance-from-ori position via
        %original angle and 'aligned' mapping of array indices
        
        %1) get absolute angle;         
        Angle=mod(atan2(ym-Chromosome.yCOM,xm-Chromosome.xCOM),2*pi);
        %1a find
        
        %2) find normalized distance vs. equi-angle
        [~,orig_idx]=min(abs(Aligned.Orig.AngleAsUsed-Angle)); %'ori' index
        SpatialPosition=Aligned.Orig.NormDist(orig_idx);  %distance;
        
        %3 find normilized, optionally corrected genomic position vs. equi-dist   
        [~,dist_idx]=min(abs(Aligned.Dist.NormAxis-SpatialPosition)); %'ori' index
        GenomicPosition=Aligned.Dist.NormCumDensityMarkercorrected(dist_idx); %genomic position
        dum=1;
        
function [Header,Headershort]=Build_cluster_Excel_Header(ClusterMaxNo);
     %this function builds headers for the table
       %Labels (case and space sensitive!):
   General = [{'index'} ,...
           {'label'}, ....
           {'cluster number'}];     
   ParamsNames =[...      
        {'content'},...     %1
        {'diameter_2Rg'},...     %2
        {'diameter_ContourBW_Equiv'},...     %2
        {'diameter_ContourBW_2Rg'},...     %2        
        {'area'},...  %3
        {'density'},...  %4
        {'xpos'},...  %5
        {'ypos'},...  %6
        {'1Ddistpos'},...  %7
        {'1DBPpos'},...  %8
        {'Psf_used'},...  %8
        ];
        
        HeaderRow1=cell(0,0); 
        HeaderRow2=cell(0,0);
        
        %1) build first part-----------------------------------------------
        GE=length(General);
        for sl=1:GE
            HeaderRow1=[HeaderRow1 {'General parameters'}];    
        end
        HeaderRow2=General;
      
        %2) build second part-----------------------------------------------
        LP=length(ParamsNames);
        for ss=1:ClusterMaxNo
            ClusterName=strcat('Cluster',num2str(ss,'%02.0f'));
            for pp=1:LP
                HeaderRow1=[HeaderRow1 cellstr(ClusterName)];
                HeaderRow2=[HeaderRow2 ParamsNames{pp}];
            end
        end
        Header=[HeaderRow1; HeaderRow2];
      dum=1;
      
      Headershort=[General(1:2) ParamsNames]; 
        
function ClusterGenprops=Get_ClusterGenprops(ClusterGenprops,ThisCellClusterTable,chro_pic,clus_pic,CellSpecs)
    %get some info on clusterfit per cell, such as off-percentage
    [clusterno,~]=size(ThisCellClusterTable);
    clus_pic=clus_pic-min(clus_pic(:));
    clus_pic=clus_pic/sum(clus_pic(:))*100;
    chro_pic=chro_pic-min(chro_pic(:));
    chro_pic=chro_pic/sum(chro_pic(:))*100;
    
    subplot(2,2,1); pcolor(chro_pic); shading flat
    subplot(2,2,2); pcolor(clus_pic); shading flat
    subplot(2,2,3); pcolor((clus_pic-chro_pic).^2); shading flat
    
    sq_diff=(sum(sum(abs(chro_pic-clus_pic))));
    %percentage uncovered by fit   
   
    idx=CellSpecs(1);
    cell_label=CellSpecs(2); 
    ClusterGenprops.cell_label(idx)=cell_label;
    ClusterGenprops.clusterno(idx)=clusterno;
    ClusterGenprops.off_perc_sq(idx)=sq_diff;
    
    dum=1;
        
%      if 0
%                 %build equivalent of avareaged spaghetticurve, but now by
%                 %cluster centers
%                 stps=25;
%                 ClusterWeightedPosAxis=(linspace(0,100,stps))';
%                 ClusterWeightedPosBP=zeros(stps,1)+NaN;
%                 ClusterWeightedPosDS=zeros(stps,1)+NaN;
%                 %build 'ClusterWeightedPosCurveBP'
%                 DS_idx=ceil((ThisCellClusterTable(:,7)+0.001)/100*(stps-1));
%                 BP_idx=ceil((ThisCellClusterTable(:,8)+0.001)/100*(stps-1));
%                 Wgt=ThisCellClusterTable(:,1);
%                 if ~isnan(sum(BP_idx)), ClusterWeightedPosBP(BP_idx)=Wgt; end
%                 if ~isnan(sum(DS_idx)), ClusterWeightedPosDS(DS_idx)=Wgt; end
% 
%                 if 0
%                     figure;
%                     stem(ClusterWeightedPosAxis,ClusterWeightedPosBP,'o-');
%                     [~]=ginput(1);
%                     close(gcf);
%                 end
%             end   
%         
        
        