function P110_cluster_replisome_I_vs_D(batchrunindex)

%analysis 1:
%-per image
%-obtain a proper replisome spot 
%-obtain cluster position and brightness;
%scatter plot intnesity of cluster vs distance to spot

close all;

if nargin<1,batchrunindex=20;end
initval=A000_Repli_Init(batchrunindex);
disp(initval.expi);
MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);

load(strcat(MatFilePath,strcat('__MovieList.mat')));
[~,LC]=size(MovieList);
int_vs_dist=[];
nearcount=0;
randomcount=0;
for ii=1:LC     
    CellName=char(MovieList(ii).CellNames); 
    disp(strcat('Done',num2str(ii),CellName,'Analysis..', num2str(LC-ii+1), 'cells to go'));
    LF=length(MovieList(ii).CellFrames);     %per frame    
    for jj=1:LF       
        %info files earlier analysis
        CellFrameName=char(MovieList(ii).CellFrames(jj)); %disp(CellFrameName);     
        CellInfoPath=strcat(MatFilePath,'ResultsOfCell',CellFrameName);
        load(strcat(CellInfoPath,'_Spots.mat'));
        load(strcat(CellInfoPath,'_Clusters.mat'));
        load(strcat(CellInfoPath,'_Cellshape.mat'));                
        [prx,pry]=get_repli_info(All_labels,'strong_ones');
        [cx,cy,cI]=get_cluster_info(Clusters,Cell);
         
         %1) pick one replisome and find nearest cluster 
         %Get_Near_clusters(prx,pry);
         LR=length(prx); 
         if~isnan(prx)
             nearrepli_idxes=[];
             for ri=1:LR  %for all replisomes
                 if ~isnan(cx)&~isempty(cx)  %if any near clusters
                     nearcount=nearcount+1;
                     %distance to clusters
                     dd_repli_clus=((cx-prx(ri)).^2+(cy-pry(ri)).^2).^0.5;  
                     [d_min,ix_nearrepli]=min(dd_repli_clus); 
                     nearrepli_idxes=[nearrepli_idxes; ix_nearrepli];
                     cI_min=cI(ix_nearrepli);          %nearest dist&I
                     cdx_min=cx(ix_nearrepli)-prx(ri);
                     cdy_min=cy(ix_nearrepli)-pry(ri);
                     
                     
                     Neigbours.repli_1stclus_drc_x(nearcount)=cdx_min;
                     Neigbours.repli_1stclus_drc_y(nearcount)=cdy_min;
                     Neigbours.repli_1stclus_distance(nearcount)=d_min;
                     Neigbours.repli_1stclus_intensity(nearcount)=cI_min;
                 end   
             end 
         end 
         %2 'random': pick one cluster and find nearest cluster
         %now exclude formerly found clusters
          LCL0=length(cx); idxes=1:LCL0;
          notformerlyfound=find(~ismember(idxes, nearrepli_idxes));
          cx_nff=cx(notformerlyfound);
          cy_nff=cy(notformerlyfound);
          cI_nff=cI(notformerlyfound);
          LCL=length(cx_nff);
         if~isnan(cx_nff)            
             for ci=1:LCL
                 if ~isnan(cx_nff)&~isempty(cx_nff)  %if any near clusters
                     
                     %distance to other clusters
                     dd_clus_clus=((cx_nff-cx_nff(ci)).^2+(cy_nff-cy_nff(ci)).^2).^0.5;
                     nonself=find(dd_clus_clus>0);
                     if~isempty(nonself)
                         randomcount=randomcount+1;
                         dd_nonself=dd_clus_clus(nonself);
                         cI_nonself=cI_nff(nonself);
                         cx_nonself=cx_nff(nonself);
                         cy_nonself=cy_nff(nonself);
                         
                         
                         [d_min,ix]=min(dd_nonself);                       
                         cI_min=cI_nonself(ix);          %nearest dist&I
                         cdx_min=cx_nonself(ix)-cx_nff(ci);
                         cdy_min=cy_nonself(ix)-cy_nff(ci);
                         Neigbours.clus_1stclus_dcc_x(randomcount)=cdx_min;
                         Neigbours.clus_1stclus_dcc_y(randomcount)=cdy_min;
                         Neigbours.clus_1stclus_distance(randomcount)=d_min;
                         Neigbours.clus_1stclus_intensity(randomcount)=cI_min;
                     end
                 end   
             end 
         end 
    end
end



plot1_x=initval.nmperpixel/1000*Neigbours.repli_1stclus_distance;
plot1_y=Neigbours.repli_1stclus_intensity;
plot2_x=initval.nmperpixel/1000*Neigbours.clus_1stclus_distance;
plot2_y=Neigbours.clus_1stclus_intensity;
plot3_x=initval.nmperpixel/1000*Neigbours.repli_1stclus_drc_x;
plot3_y=initval.nmperpixel/1000*Neigbours.repli_1stclus_drc_y;
plot3_mrk=round(Neigbours.repli_1stclus_intensity/10+3);

plot4_x=initval.nmperpixel/1000*Neigbours.clus_1stclus_dcc_x;
plot4_y=initval.nmperpixel/1000*Neigbours.clus_1stclus_dcc_y;
plot4_mrk=round(Neigbours.clus_1stclus_intensity/10+3);

binax=linspace(1,100,50);
hist1=hist(plot1_y,binax); hist1=hist1/sum(hist1)*100;
hist2=hist(plot2_y,binax); hist2=hist2/sum(hist2)*100;

subplot(2,2,1);
    plot(plot1_x,plot1_y,'ro','MarkerSize',3,'MarkerFaceColor','r'); hold on;
    plot(plot2_x,plot2_y,'bo','MarkerSize',3); hold on;   
    title([initval.expi,':cluster intensity vs replisome distance']);
    xlabel('distance, microns');
    ylabel('content, %of image');
    xlim([0 1.5]); ylim([0 100]);
    legend('1st cluster next to replisome','other clusters');
    
subplot(2,2,3);
    L3=length(plot3_x);
    for i3=1:L3
        plot(plot3_x(i3),plot3_y(i3),'ro','MarkerSize',plot3_mrk(i3)); hold on;
    end
    L4=length(plot4_x);
    for i4=1:L4
        plot(plot4_x(i4),plot4_y(i4),'bo','MarkerSize',plot4_mrk(i4)); hold on;  
    end
    title([initval.expi,':cluster intensity vs replisome distance']);
    axis square
    xlabel('distance, microns');
    ylabel('distance, microns');
     xlim([-1 1]); ylim([-1 1]);
    legend('replisome-1st cluster','cluster-1st cluster');
subplot(2,2,2);
    bar(binax,hist1,'r');
    xlabel('content, %of image');
    ylabel('occurrence, %');
    xlim([0 100]); %ylim([0 100]);
    legend('1st cluster next to replisome');
subplot(2,2,4);
    bar(binax,hist2,'b');
    xlabel('content, %of image');
    ylabel('occurrence, %');
    xlim([0 100]); %ylim([0 100]);
    legend('other clusters');    
    
if 0
subplot(1,2,2);
        Ydata=plot1_y;
        Xdata=plot1_x;
        lobinY=1;    %nm
        hibinY=100;   
        binsY=100;
        binax_Y=linspace(lobinY,hibinY,binsY);

        lobinX=0;     %seconds
        hibinX=10;
        binsX=100;
        binax_X=linspace(lobinX,hibinX,binsX);
        [YY,XX]=meshgrid(binax_Y, binax_X);
        ExpandData=Expand_Data(Ydata,Xdata,lobinY,hibinY,lobinX,hibinX);
        YX_hist=hist3([ExpandData],{binax_Y binax_X});
        YX_hist=YX_hist/max(YX_hist(:));
        [CC,h] =contour(XX,YY,YX_hist',25); colormap hot; shading flat; hold on;        
        xlabel('distance, microns');
        ylabel('content, %of image');
        xlim([1 10]); ylim([1 100]);
end
 pause(0.01);           
saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_P110_cluster_replisome_distances.jpg'),'jpg');
save(strcat(initval.resultpath,initval.DirSep,initval.expi,'_P110_cluster_replisome_distances.mat'),'Neigbours');

pause(0.01); 


function [clusters_x,clusters_y,clusters_I]=get_cluster_info(Clusters,Cell)
        %Clusters:
        % 'spotprops'
        % 'COM_X'
        % 'COM_Y'
        % 'C_perc'
        % 'EquivDiameter'
        % 'Area'
        % 'Density'
        % 'psf_used'
        % 'clustermask'
        [~,LC]=size(Clusters);
        if LC>0
        for ci=1:LC   %get clusters
            clusters_x(ci)=Clusters(ci).COM_X;
            clusters_y(ci)=Clusters(ci).COM_Y;
            clusters_I(ci)=Clusters(ci).C_perc;
        end
        else
            clusters_x=NaN;
            clusters_y=NaN;
            clusters_I=NaN;
        end
        clusters_I=100*clusters_I; %*Cell.N_est;
        
   function [repli_x,repli_y]=get_repli_info(All_labels,seltype);
        %stored replisome info                   
        %All_labels.Rfp:
        % 'OK_nicepic'
        % 'OK_pair'
        % 'OK_pair_dist'
        % 'OK_pair_idx'
        % 'OK_strong'
        % 'spotContent'
        % 'spotX'
        % 'spotY'
        % 'spot_ori_index'
        % 'frameindex'
        badrepli=0;
        if isfield(All_labels.Rfp, 'OK_pair')
            switch seltype
                case 'strong_ones', sel=find(All_labels.Rfp.OK_strong==1);  %proper spots
                case 'paired', sel=find(All_labels.Rfp.OK_pair==1);  %proper spots   
            end
            if ~isempty(sel)
            repli_x=All_labels.Rfp.spotX(sel);
            repli_y=All_labels.Rfp.spotY(sel);
            badrepli=0;
            else
               badrepli=1; 
            end
        else badrepli=1;
        end
        if badrepli==1;
            repli_x=NaN;
            repli_y=NaN;
        end
        
    function  ExpandData=Expand_Data(Y, X,lobinY,hibinY,lobinX,hibinX)
    %This function renders sparse data points clouds usable for contour histogram use;
    %it does so by creating gaussian distributed point clouds around every
    %data point. Effectively, this smooths a contour hisogram which would
    %otherwise be very fragmented. The smoothing width covers a fifth of
    %the standard deviation of step sizes or dwell times. JacobKers 2017
    LD=length(X);
    ExpandFactor=500;
    loX=min(Y); hiX=max(Y); rangeX=nanstd(Y);
    loT=min(X); hiT=max(X); rangeT=nanstd(X);
    ExpandRangeX=rangeX/5;
    ExpandRangeT=rangeT/5;
    
    ExpandX=[]; ExpandT=[];
    for ii=1:LD  
        scatX=Y(ii)+ExpandRangeX*randn(ExpandFactor,1);
        scatT=X(ii)+ExpandRangeT*randn(ExpandFactor,1);
        ExpandX=[ExpandX; scatX];
        ExpandT=[ExpandT; scatT];
    end
    
    %cropping
    sel1=find((ExpandT<lobinX)|ExpandT>hibinX); ExpandT(sel1)=NaN;
    sel2=find((ExpandX<lobinY)|ExpandX>hibinY); ExpandX(sel2)=NaN;
    ExpandData=[ExpandX ExpandT];
    
        