function [N_clust,AllClusters,Clusters,GenProps]=TCl_201809_get_clusters(im,initval)
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

%close all;
savit=1;
shoit=1;
if nargin<2
    initval.Psf_est=2.4;
    initval.skipsmallspots=2/100;  %max fraction of spot to accept
    initval.Separation=2.5;
    initval.contourtreshold=0.15;
    initval.outpath=pwd;
    %margin to re-connect spots, in units of psf.
    %Larger than two related to non-equal sized spots
    hand_load=1;
    if~hand_load
        [CellName,im]=TCl_201809_build_cluster_testpattern;  
    else
        close all;
        curp=pwd;
        cd('C:\Users\jkerssemakers\Dropbox\CD_recent\BN_CD16_Sandro\Paper2019\Replisome_Cluster-Colocalization');
        [FileName,PathName] = uigetfile('.tif');
        CellName=FileName(1:end-4);
        im=double(imread(FileName))
        im=im-min(im(:));
        cd(curp);
        figure(88);
            plot_pic=-im'+max(im(:)); 
            pcolor(plot_pic); colormap bone; shading flat; 
            axis equal; axis off; hold on;
    end
    initval.outpath=pwd;
    initval.CellName=CellName;
    imoutdir=initval.outpath;
    MatFileDir=initval.outpath;
else
    imoutdir=strcat(initval.outpath,'\Results_clus\');
    MatFileDir=strcat(initval.outpath,'\Results_clus\');
end
Psf_meas=initval.Psf_est;


AllCellsAllClusters=[];
AllClusters=[];
ClusterMaxNo=0;
im_bc=TCl_201809_remove_cluster_background(im,'Min');
    AllSpotProps=TCl_201809_peel_spots_from_image(im_bc,Psf_meas,1,initval.skipsmallspots, 0); 
    AllSpotProps=TCl_201809_clean_spot_collection(AllSpotProps,im_bc,Psf_meas);      
    [rr,cc]=size(AllSpotProps);        
    if rr>1  %if found spots
        ClusterBasics=Find_Clusters(AllSpotProps, initval);
        [~,CL]=size(ClusterBasics);
        %Build single cluster contours
        [Clusters,ClusterPropsRow,ThisCellClusterTable,GenProps]=GetClusterDetails(ClusterBasics,im_bc,initval);
        
        [~,clusterno]=size(Clusters);
           
        ClusterMaxNo=max([ClusterMaxNo clusterno]); %needed for table header

        %add to full per-cell-table, elongate if necessary
        [rr,cc]=size(AllCellsAllClusters);
        newrowentry=[clusterno ClusterPropsRow]; le=length(newrowentry);
        cc_new=max([cc,le]);        
        newtable=zeros(rr+1,cc_new); 
        newtable(1:rr,1:cc)=AllCellsAllClusters;
        newtable(rr+1,1:le)=newrowentry;
        AllCellsAllClusters=newtable;

        %add to cluster-sorted table
        AllClusters=[AllClusters; ...
                    [ThisCellClusterTable]];
        
        N_clust=CL;
        
        
%% PLOT MENU---------------------------------------------------------------
%% plot the cromosome pattern
    %things we plot per cell
    if savit
%          subplot(1,2,1);
%          plot(GenProps.AllContoursY,GenProps.AllContoursX,'y-','LineWidth',2); hold on; 
%          subplot(1,2,2);
         plot(GenProps.AllContoursY,GenProps.AllContoursX,'y-','LineWidth',2); hold on; 
        addspots=0;
        if addspots
            figure(88);
            plot_pic=-im_bc'+max(im_bc(:)); 
            pcolor(plot_pic); colormap bone; shading flat; 
            axis equal; axis off; hold on;
            PlotProps=AllSpotProps;   
            perc=num2str(round(100*PlotProps(end,7)));   
            scale=24/max(PlotProps(:,6));   
            plot(GenProps.AllContoursY,GenProps.AllContoursX,'y-','LineWidth',2); hold on;
            scale=24/max(PlotProps(:,6));
            [LS,~]=size(PlotProps);
            for sp=1:LS
                mrksize=ceil(PlotProps(sp,6)*scale);
                x=PlotProps(sp,3);
                y=PlotProps(sp,4);

                plot(y,x,'wo','MarkerSize',max([mrksize 1]),'LineWidth',2);
                hold on; 
                %text(y+1,x+1,num2str(sp),'BackgroundColor','w');
                %axis equal;
            end
        end
        title(initval.CellName)
        saveas(gcf,strcat(imoutdir,'clusters_', initval.CellName,'.jpg')); 
        pause(0.1);
    end
 end
    
%% wrapping up for saving: cluster report
 if savit
     save(strcat(MatFileDir,initval.CellName,'_Clusters.mat'),'Clusters','im');  
     
     if 0
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
    erasersheet=NaN*zeros(1E3,2*cc);
    xlswrite(strcat(MatFileDir,initval.CellName,'_A050_Cell_ClusterReport.xlsx'),Header,'Sheet1','A1');
    xlswrite(strcat(MatFileDir,initval.CellName,'_A050_Cell_ClusterReport.xlsx'),erasersheet,'Sheet1','A3');
    xlswrite(strcat(MatFileDir,initval.CellName,'_A050_Cell_ClusterReport.xlsx'),AllCellsAllClusters,'Sheet1','A3');
    save(strcat(MatFileDir,initval.CellName,'_A050_Cell_ClusterReport.mat'),'AllCellsAllClusters');

    [rr2,cc2]=size(AllClusters);
    erasersheet2=NaN*zeros(1E4,2*cc2);

    xlswrite(strcat(MatFileDir,initval.CellName,'_A050_ClusterColumnReport.xlsx'),Headershort,'Sheet1','A1');
    xlswrite(strcat(MatFileDir,initval.CellName,'_A050_ClusterColumnReport.xlsx'),erasersheet2,'Sheet1','A2');
    xlswrite(strcat(MatFileDir,initval.CellName,'_A050_ClusterColumnReport.xlsx'),AllClusters,'Sheet1','A2');
    end
    
    save(strcat(MatFileDir,initval.CellName,'_A050_ClusterColumnReport.mat'),'AllClusters');
 end

function [pairsX,pairsY]=Find_Pairs(sp,PlotProps,initval);
%This function finds neighbouring points of a specified centre point within a
%specified range. It outputs the connections between the centre points and
%these neighbours as coordinate pairs. JacobKers2016
    x0=PlotProps(sp,3);
    y0=PlotProps(sp,4);
    
    sel=find(PlotProps(:,1)~=sp);
    otherspots=PlotProps(sel,:);
    otherx=otherspots(:,3);
    othery=otherspots(:,4);
    rr=((otherx-x0).^2+(othery-y0).^2).^0.5;

    near_ones=find(rr<initval.Separation*initval.Psf_est);
    pairsX=[otherx(near_ones)' ; 0*otherx(near_ones)'+x0];
    pairsY=[othery(near_ones)' ; 0*othery(near_ones)'+y0];

function Clusters=Find_Clusters(SpotProps, initval);
% 'Use this section for a Quicksheet'
%------------------------------------------------------------
    % sort spots by brightness; pick brightest one
    % find near ones for brightmost one from leftover list
    % for the ones found, keep finding new near ones until no change; 
   
        % repeat for leftovers until no leftovers
        %JacobKers 2017 ----------------------------------------
        %AllSpotProps=[
        % 1 spotcount 
        % 2 Peak 
        % 3 Xpos 
        % 4 Ypos 
        % 5 Psf 
        %6 ThisSpotFraction 
        %7CoveredFraction 
        %8 RelChange]];
%------------------------------------------------------------[JK16]
% 'End of Quicksheet section'
 psf=initval.Psf_est;
 left_ones=SpotProps;
 Clusters=struct([]);
 N_clust=0;
 ClusterTable=[];
 while ~isempty(left_ones);    
    %build a new cluster, starting with the brightest leftover spot
    N_clust=N_clust+1;
    Bri=left_ones(:,2); [Brisort,idx]=sort(Bri,'descend');
    left_ones=left_ones(idx,:);
    thiscluster=left_ones(1,:);
    %start spot
     left_ones=left_ones(2:end,:);  %pick from stock
    xc0=thiscluster(3); yc0=thiscluster(4);
    
    cluster_grows=1;
    while cluster_grows
        [LC,~]=size(thiscluster);
        new_ones=[];
        for ii=1:LC  
            %for all existing elements of this cluster, 
            %,find near ones in the leftover stock
            x0=thiscluster(ii,3);
            y0=thiscluster(ii,4);
            a0=thiscluster(ii,6); %fraction of whole pattern
            
            leftx=left_ones(:,3);
            lefty=left_ones(:,4);
            lefta=left_ones(:,6); %fraction of whole pattern
         
            ispartofcluster=Get_Separation_conditions(x0,y0,a0,leftx,lefty,lefta,initval);
                  
            near_ones_idx=find(ispartofcluster); %overlapping
            left_ones_idx=find(~ispartofcluster); 
          
            
            if ~isempty(near_ones_idx)
                new_ones=[new_ones; left_ones(near_ones_idx,:)]; %add
                left_ones=left_ones(left_ones_idx,:);  %shrink remaining                 
            end
        end
        if ~isempty(new_ones) 
            thiscluster=[thiscluster; new_ones]; %add new ones to cluster
            
            if 0
            %subplot(1,2,2);
            plot(thiscluster(:,4),thiscluster(:,3),'ko'); hold on;
            plot(yc0,xc0,'ro');
            plot(new_ones(:,4),new_ones(:,3),'bo'); hold on;
            axis equal;
             axis([0 70 0 70]);
            [~]=ginput(1);
            hold off;
            end
            
        else
            cluster_grows=0;  %...or stop this cluster
        end
    end
    %N_clust=N_clust+1;
    %add values to cluster struct; go to next
    Clusters(N_clust).spotprops=thiscluster;    
 end

 function ispartofcluster=Get_Separation_conditions(x0,y0,a0,leftx,lefty,lefta,initval);
            %define an amplitude-dependent separation condition. The weaker
            %a spot, the more it is allowed to connect. As a definition for
            %connection, we take when the minimum between two unequal
            %1D Gaussians dissapears, estimated to be where the larger
            %gaussian equals half of the peak height of the smaller one.         
            
            %connect if: 
            %-close enough
            %-correct for differences
            %-both strong enough (avoid bridging valleys)
            psf=initval.Psf_est;
            LL=length(lefta);
            rr=((leftx-x0).^2+(lefty-y0).^2).^0.5;            
           
            dif_cor=(abs(2*log(lefta./(2*a0)))).^0.5;
            
            maxsep=initval.Separation*psf*dif_cor;                        
            maxsep=(min([maxsep'; (ones(1,LL))*2*initval.Separation*psf]))'; %clip  
            maxsep=initval.Separation*psf; 
            aretheothersclose=rr<maxsep;
            
            isthisstrongspot=a0>2*initval.skipsmallspots;
            isotherstrongspot=lefta>2*initval.skipsmallspots;
            
            ispartofcluster=isthisstrongspot&aretheothersclose;

function [Clusters,ThisCellClusterRow,ThisCellClusterTable,GenProps]=GetClusterDetails(Clusters,orim,initval)
% 'Use this section for a Quicksheet'
%------------------------------------------------------------
    % this function works trhough all clusters in a cell:
    % builds a single-cluster image and obtains properties for this.
     %input: AllSpotProps, collection of contributing Gaussian spots
        % 1 spotcount 
        % 2 Peak 
        % 3 Xpos 
        % 4 Ypos 
        % 5 Psf 
        %6 ThisSpotFraction 
        %7CoveredFraction 
        %8 RelChange]];
    % It also builds a clustercontour(level) using the FWHM lelvel 
    % points (or other levels)
        %Clusters
        %Contours
        %ThisCellClusterTable contains sorted  'CL' rows of:
        %order content% area radius density xpos ypos 1Ddistpos 1DBPpos
%------------------------------------------------------------[JK16]
% 'End of Quicksheet section'
 
    [rr,cc]=size(orim);
    [XX,YY]=meshgrid(1:cc,1:rr);    
    [~,CL]=size(Clusters);
    AllContoursX=zeros(50,CL);
    AllContoursY=zeros(50,CL);
    ThisCellClusterTable=zeros(CL,9);
    %contains 'CL' rows of:
    %order content% area radius density xpos ypos 1Ddistpos 1DBPpos
    ReconstructIm=0*orim;   
     for ii=1:CL  %FOR ALL CLUSTERS
             
        ThisClusterSpots=Clusters(ii).spotprops;  %components of cluster
        thisclusterim=GetSingleClusterImage(ThisClusterSpots,orim); 
        ReconstructIm=ReconstructIm+thisclusterim;
        C_perc=sum(ThisClusterSpots(:,6));
        psf_used=ThisClusterSpots(1,5);
        
        %cluster parameters  
        [x_com,y_com,~,~,rad_gyr]=JKD2_IM_calculate2Dmoment_extended(thisclusterim);
        if 0 %strcmp(cellno,'200020')
            pcolor(thisclusterim'); shading flat, colormap bone;
            title(['Cluster' num2str(ii)]);
             axis equal; axis tight; axis xy; hold on;
            [~]=ginput(1);
        end   
        C_tres=initval.contourtreshold;
        [contourX,contourY,ClusterShapeProps,clustermask]=GetClusterShapeProps(thisclusterim,orim,x_com,y_com,C_tres); 
        
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
                          psf_used...
                          ];
            %order content% area radius density xpos ypos 1Ddistpos 1DBPpos
        ThisCellClusterTable(ii,:)=thisclusterprops;
      
        %add properties to 'Cluster' structure
        Clusters(ii).COM_X=x_com;  %[x,y];
        Clusters(ii).COM_Y=y_com;  %[x,y]; 
        Clusters(ii).C_perc=C_perc;                       
        Clusters(ii).EquivDiameter=ClusterShapeProps.EquivDiameter;
        Clusters(ii).Area=ClusterShapeProps.Area; 
        Clusters(ii).Density=ClusterShapeProps.Density;
        Clusters(ii).psf_used=psf_used; 
        Clusters(ii).clustermask=clustermask;    
     end
        GenProps.orim=orim;
        GenProps.clus_im=ReconstructIm;
        GenProps.AllContoursX=AllContoursX;
        GenProps.AllContoursY=AllContoursY;
     
     %Sort cluster table by brightest cluster first, transform  to one row
     %content% radius area density xpos ypos 1Ddistpos 1DBPpos
     [rr,cc]=size(ThisCellClusterTable);
     [~,idx]=sort(ThisCellClusterTable(:,1),'descend');
     ThisCellClusterTable=ThisCellClusterTable(idx,:);
     ThisCellClusterRow=reshape(ThisCellClusterTable',1,rr*cc); %repeat=8;
      
  function [contourX,contourY,ClusterShapeProps,BW]=GetClusterShapeProps(thisclusterim,orim,xm,ym,cutoff);
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
                [contourX,contourY]=TCl_201809_B002_EqualizeAlongContour(contourX,contourY,49); 
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
      Headershort=[General(1:2) ParamsNames]; 
        

    
   
      
        
        
        
        
        