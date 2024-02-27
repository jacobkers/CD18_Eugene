function A055_WF_2DClusterAnalysis_Standalone(initval)
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
%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------

close all;
action.runanalysis=1;
action.runplots=1;

if nargin<1
    usr='Jacob', batchrunindex=0;
    initval=A000_Repli_Init(batchrunindex,usr);
end

Psf_meas=initval.Psf_est;
allframes=length(initval.Cell_Labels);

%% set up output directory
imoutdir=strcat(initval.pth_repli,'CellImages_Clusteranalysis',initval.DirSep);  
if isdir(imoutdir), rmdir(imoutdir,'s');  end
mkdir(imoutdir);
MatFileDir=strcat(initval.pth_repli,'ResultsPerCellMatlab',initval.DirSep); 


%% loop 1: cluster analysis
AllCellsAllClusters=[];
AllClusters=[];
goodcount=0;
ClusterMaxNo=0;
if action.runanalysis
    for jj=1:allframes 
        %naming and labelling stuff
        goodcount=goodcount+1; 
        cellno=char(initval.Cell_Labels{jj}); 
        NumCellLabel=BuildNumericCellLabel(cellno);
        CellSpecs=[goodcount NumCellLabel]; 
        CellName=strcat('ResultsOfCell',cellno); 
        disp(strcat('Program:A055_experiment:',initval.expi,':',CellName,'ClusterAnalysis..', num2str(allframes-jj+1), 'cells to go'));    
          
        channel_stack=get_channel_stack_from_cropcode(cellno,initval);
        [~, cellmask, chro_pic, ~, ~]=get_channel_ID(channel_stack, initval);

         %bit of refinement
        chro_pic=(chro_pic-min(chro_pic(:))).*cellmask;
        celledge_pic = bwmorph(cellmask,'remove');     
        
        AllSpotProps=PeelblobsFromImage(chro_pic,Psf_meas,0.98,0.03, 0,initval.cluster_sep_sigs);
        %AllSpotProps=F006_CleanSpots(AllSpotProps,chro_pic,Psf_meas);      
        
        %next, group these into clusters
        [rr,~]=size(AllSpotProps);        
        if rr>1
        ClusterBasics=Find_Clusters(AllSpotProps, initval);
        [Clusters,ClusterPropsRow,ThisCellClusterTable,AllContours,ReconstructIm]=...
            GetClusterDetails(ClusterBasics,chro_pic);
        [~,clusterno]=size(Clusters);
        
        %add to full per-cell-table, elongate if necessary
        ClusterMaxNo=max([ClusterMaxNo clusterno]); %needed for table header     
        [rr,cc]=size(AllCellsAllClusters);
        newrowentry=[CellSpecs clusterno ClusterPropsRow]; le=length(newrowentry);
        cc_new=max([cc,le]);        
        newtable=zeros(rr+1,cc_new); 
        newtable(1:rr,1:cc)=AllCellsAllClusters;
        newtable(rr+1,1:le)=newrowentry;
        AllCellsAllClusters=newtable;

        %add to cluster-sorted table
        AllClusters=[AllClusters; ...
                    [repmat(CellSpecs,clusterno,1)...
                     ThisCellClusterTable]];                                                                          
        else
        Clusters=struct([]);        
        end
        save(strcat(MatFileDir,CellName,'_Clusters.mat'),'Clusters','AllContours','chro_pic','celledge_pic'); 
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
        disp('Saving Excel.....');
        xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A050_Cell_ClusterReport.xlsx'),Header,'Sheet1','A1');
        xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A050_Cell_ClusterReport.xlsx'),erasersheet,'Sheet1','A3');
        xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A050_Cell_ClusterReport.xlsx'),AllCellsAllClusters,'Sheet1','A3');
        save(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A050_Cell_ClusterReport.mat'),'AllCellsAllClusters');

        [rr2,cc2]=size(AllClusters);
        erasersheet2=NaN*zeros(10E4,2*cc2);
        xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A055_ClusterColumnReport.xlsx'),Headershort,'Sheet1','A1');
        xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A055_ClusterColumnReport.xlsx'),erasersheet2,'Sheet1','A2');
        xlswrite(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A055_ClusterColumnReport.xlsx'),AllClusters,'Sheet1','A2');
        save(strcat(initval.pth_repli,initval.DirSep,initval.expi,'_A055_ClusterColumnReport.mat'),'AllClusters');

end

%% plot loop
if action.runplots
    for jj=1:allframes
        if initval.showplots
            set(figure(1), 'visible','on')
        else
            set(figure(1), 'visible','off')
        end
        cellno=char(initval.Cell_Labels{jj}); 
        CellName=strcat('ResultsOfCell',cellno); 
        
        channel_stack=get_channel_stack_from_cropcode(cellno,initval);
        [~, cellmask, chro_pic, ~, ~]=get_channel_ID(channel_stack, initval);

         %bit of refinement
        chro_pic=(chro_pic-min(chro_pic(:))).*cellmask;
        celledge_pic = bwmorph(cellmask,'remove');     
        
        disp(strcat(CellName,'Plotting..', num2str(allframes-jj+1), 'cells to go'));
        load(strcat(MatFileDir,CellName,'_Clusters.mat'),'Clusters','AllContours','chro_pic','celledge_pic');          
            whereedgecell=find(celledge_pic~=0);
              
            load(strcat(MatFileDir,CellName,'_Spots.mat'));  
            
            subplot(2,2,1);                          
            AllContoursX=AllContours.X;
            AllContoursY=AllContours.Y;
                                        
            if initval.analyze_ori_ter
                Cfp=All_labels.Cfp;
                Cfp_pic=All_labels.Cfp_pic;
                terX=round(Cfp.spotX); terY=round(Cfp.spotY);
                plot(terX+0.5,terY+0.5, 'bo','MarkerSize',10, 'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;
                plot(terX+0.5,terY+0.5, 'bo','MarkerSize',10, 'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;
            end
                if initval.analyze_replisome | initval.analyze_ori_ter
                Rfp=All_labels.Rfp;  
                Rfp_plotpic=All_labels.Rfp_pic;
                
                
                
                pcolor(Rfp_plotpic); shading flat, colormap bone; axis equal; axis off; hold on;   
                pause(0.02);
                [cellcontourX,cellcontourY]=get_smooth_contour(celledge_pic,100);
                %Rfp_plotpic(whereedgecell)=max(Rfp_plotpic(:)); 
                plot(cellcontourX,cellcontourY,'w-','Linewidth',1); hold on;
                wk=find(Rfp.OK_strong==0);  %weak spots
                if length(wk)>0,
                    plot(Rfp.spotX(wk),Rfp.spotY(wk), 'yx', 'MarkerSize', 6); 
                end
                strng=find(Rfp.OK_strong==1); %strongspots
                plot(Rfp.spotX(strng),Rfp.spotY(strng), 'wo', 'MarkerSize', 8); 
                %build pairs
                sel=find(Rfp.OK_pair==1);
                LP=length(sel); pairX=zeros(LP,2); pairY=zeros(LP,2);
                for pp=1:LP
                    id1=sel(pp); 
                    id2=Rfp.OK_pair_idx(id1);
                    pairX(pp,:)=[Rfp.spotX(id1) Rfp.spotX(id2)];
                    pairY(pp,:)=[Rfp.spotY(id1) Rfp.spotY(id2)];           
                end       
                plot(pairX',pairY','wo-', 'MarkerSize',4,'MarkerFaceColor','w');
                title('rfp');
                hold off
           end
       subplot(2,2,2);
            chroplot3_pic=chro_pic;
            %chroplot3_pic(whereedgecell)=max(chroplot3_pic(:)); 
            pcolor(chroplot3_pic); shading flat, colormap bone; axis equal; hold on;
            pause(0.02);
             [cellcontourX,cellcontourY]=get_smooth_contour(celledge_pic,100);
                %Rfp_plotpic(whereedgecell)=max(Rfp_plotpic(:)); 
                plot(cellcontourX,cellcontourY,'w-','Linewidth',1); hold on;
            axis off; hold off;   
            title('chromosome');
      subplot(2,2,3);
            if initval.analyze_replisome | initval.analyze_ori_ter
                Rfp=All_labels.Rfp; 
                chroplot_pic=-chro_pic+max(chro_pic(:));  
                %chroplot_pic(whereedgecell)=min(chroplot_pic(:)); 
                pcolor(chroplot_pic); shading flat, colormap bone; axis equal; hold on;            
                
                pause(0.02);
                [cellcontourX,cellcontourY]=get_smooth_contour(celledge_pic,100); hold on;
                %Rfp_plotpic(whereedgecell)=max(Rfp_plotpic(:)); 
                plot(cellcontourX,cellcontourY,'k-','Linewidth',1);
                contour(chro_pic,6,'k','Linewidth',1); hold on;
                axis equal; hold on; 
                load(strcat(MatFileDir,CellName,'_Spots.mat'));            
                if initval.analyze_ori_ter
                    Cfp=All_labels.Cfp;
                    Cfp_pic=All_labels.Cfp_pic;
                    terX=round(Cfp.spotX); terY=round(Cfp.spotY);
                    plot(terX+0.5,terY+0.5, 'bo','MarkerSize',10, 'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;
                    plot(terX+0.5,terY+0.5, 'bo','MarkerSize',10, 'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;

                end
                if initval.analyze_replisome | initval.analyze_ori_ter              
                    Rfp=All_labels.Rfp;  
                    wk=find(Rfp.OK_strong==0);  %weak spots
                    if length(wk)>0,
                        plot(Rfp.spotX(wk),Rfp.spotY(wk), 'yx', 'MarkerSize', 6); 
                    end
                    strng=find(Rfp.OK_strong==1); %strongspots
                    plot(Rfp.spotX(strng),Rfp.spotY(strng), 'ro', 'MarkerSize', 8); 
                    %build pairs
                    sel=find(Rfp.OK_pair==1);
                    LP=length(sel); pairX=zeros(LP,2); pairY=zeros(LP,2);
                    for pp=1:LP
                        id1=sel(pp); 
                        id2=Rfp.OK_pair_idx(id1);
                        pairX(pp,:)=[Rfp.spotX(id1) Rfp.spotX(id2)];
                        pairY(pp,:)=[Rfp.spotY(id1) Rfp.spotY(id2)];           
                    end       
                    plot(pairX',pairY','ro-', 'MarkerSize',5,'MarkerFaceColor','r');  
                end  
                plot(AllContoursX,AllContoursY,'y-','LineWidth',1); hold off;
                axis equal;
                title(strcat('Cell', num2str(cellno,'% 3.0f')));
            end
            saveas(gcf,strcat(imoutdir,'CellClusters', num2str(cellno,'% 3.0f'),'.jpg')); 
            close(gcf);
            if initval.showplots
                
                set(figure(2), 'visible','on')
            else
                set(figure(2), 'visible','off')
            end
               if initval.analyze_replisome | initval.analyze_ori_ter
                Rfp=All_labels.Rfp; 
                chroplot_pic=-chro_pic+max(chro_pic(:));  
                %chroplot_pic(whereedgecell)=min(chroplot_pic(:)); 
                pcolor(chroplot_pic); shading flat, colormap bone; axis equal; hold on;            
                pause(0.02);
                [cellcontourX,cellcontourY]=get_smooth_contour(celledge_pic,100); hold on;
                %Rfp_plotpic(whereedgecell)=max(Rfp_plotpic(:)); 
                plot(cellcontourX,cellcontourY,'k-','Linewidth',1);
                
                contour(chro_pic,6,'k','Linewidth',1); hold on;
                axis equal; hold on; 
                load(strcat(MatFileDir,CellName,'_Spots.mat'));            
                if initval.analyze_ori_ter
                    Cfp=All_labels.Cfp;
                    Cfp_pic=All_labels.Cfp_pic;
                    terX=round(Cfp.spotX); terY=round(Cfp.spotY);
                    plot(terX+0.5,terY+0.5, 'bo','MarkerSize',10, 'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;
                    plot(terX+0.5,terY+0.5, 'bo','MarkerSize',10, 'MarkerFaceColor','b','MarkerEdgeColor','k'); hold on;

                end
                if initval.analyze_replisome | initval.analyze_ori_ter              
                    Rfp=All_labels.Rfp;  
                    wk=find(Rfp.OK_strong==0);  %weak spots
                    if length(wk)>0,
                        plot(Rfp.spotX(wk),Rfp.spotY(wk), 'yx', 'MarkerSize', 6); 
                    end
                    strng=find(Rfp.OK_strong==1); %strongspots
                    plot(Rfp.spotX(strng),Rfp.spotY(strng), 'ro', 'MarkerSize', 8); 
                    %build pairs
                    sel=find(Rfp.OK_pair==1);
                    LP=length(sel); pairX=zeros(LP,2); pairY=zeros(LP,2);
                    for pp=1:LP
                        id1=sel(pp); 
                        id2=Rfp.OK_pair_idx(id1);
                        pairX(pp,:)=[Rfp.spotX(id1) Rfp.spotX(id2)];
                        pairY(pp,:)=[Rfp.spotY(id1) Rfp.spotY(id2)];           
                    end       
                    plot(pairX',pairY','ro-', 'MarkerSize',5,'MarkerFaceColor','r');  
                end  
                plot(AllContoursX,AllContoursY,'y-','LineWidth',2); hold off;
                axis equal;
                title(strcat('Cell', num2str(cellno,'% 3.0f')));
            end
            saveas(gcf,strcat(imoutdir,'CellClusters', num2str(cellno,'% 3.0f'),'_B.jpg'));
            close(gcf);
    end
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
    near_ones=find(rr<2.0*initval.Psf_est);
    pairsX=[otherx(near_ones)' ; 0*otherx(near_ones)'+x0];
    pairsY=[othery(near_ones)' ; 0*othery(near_ones)'+y0];

function Clusters=Find_Clusters(SpotProps, initval);
% 'Use this section for a Quicksheet'
%------------------------------------------------------------
    % sort spots by brightness; pick brightest one
    % find near ones for brightmost one from leftover list
    % for the ones found, keep finding new near ones until nochange; 
   
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
SeparateWidth=initval.cluster_sep_sigs;
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


function [Clusters,ThisCellClusterRow,ThisCellClusterTable,AllContours,ReconstructIm]=GetClusterDetails(Clusters,orim)
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
        %[x_com,y_com,~,~]=TrackXY_by_COM_2Dmoment(thisclusterim);
        
        
        [x_com,y_com,~,~,rad_gyr]=JKD2_IM_calculate2Dmoment_extended(thisclusterim);
        if 0 %strcmp(cellno,'200020')
            pcolor(thisclusterim'); shading flat, colormap bone;
            title(['Cluster' num2str(ii)]);
             axis equal; axis tight; axis xy; hold on;
            [~]=ginput(1);
        end                               
        [contourX,contourY,ClusterShapeProps,clustermask]=GetClusterShapeProps(thisclusterim,orim,x_com,y_com,0.2); 
        
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
     AllContours.X=AllContoursX;
     AllContours.Y=AllContoursY;
     
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
                blowup=5; 
                [rr,cc]=size(thisclusterim);
                [XX,YY]=meshgrid(1:cc,1:rr);
                if blowup>1
                    colax=linspace(1,cc,blowup*cc);
                    rowax=linspace(1,rr,blowup*rr);
                    [XXip, YYip]=meshgrid(colax,rowax);
                    thisclusterim_blw=interp2(XX,YY,thisclusterim,XXip, YYip);
                    %thisclusterim=thisclusterim_blw;
                    thisclusterim=JKD2_IM_smoothJK(thisclusterim_blw,blowup/2);
                    [rr,cc]=size(thisclusterim);
                    [XX,YY]=meshgrid(1:cc,1:rr);
                    xm=xm*blowup;
                    ym=ym*blowup;
                end
                
                BW=0*thisclusterim;
                sel=find(thisclusterim>cutoff*max(orim(:)));
                if length(sel)>1
                BW(sel)=1;   
                ClusterShapeProps = regionprops(BW, 'Area','EquivDiameter');
                ClusterShapeProps.Density=sum(thisclusterim(sel))/length(sel);
                
                [xBW,yBW,~,~,RgBW]=JKD2_IM_calculate2Dmoment_extended(1.0*BW);
    
                
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
                    xBW=NaN; 
                    yBW=NaN; 
                    RgBW=NaN; 
                end
                
                %correct for blowup
                ClusterShapeProps.xBW=xBW/blowup; 
                ClusterShapeProps.yBW=yBW/blowup; 
                ClusterShapeProps.RgBW=RgBW/blowup;
                ClusterShapeProps.Area=ClusterShapeProps.Area/blowup.^2;
                ClusterShapeProps.EquivDiameter=ClusterShapeProps.EquivDiameter/blowup;
                ClusterShapeProps.xBW=ClusterShapeProps.xBW/blowup; 
                ClusterShapeProps.yBW=ClusterShapeProps.yBW/blowup;                
                contourX=contourX/blowup;
                contourY=contourY/blowup;
                dum=1;
                
                
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
      Headershort=[General(1:2) ParamsNames]; 
        
    
      
        
        
        
        
        