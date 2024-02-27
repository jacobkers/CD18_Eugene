function TCl_201809_process_noiseplotter
% Plotdata:
    % hist_binax: [1×51 double]
    % hist_meas: [1×51 double]
    % hist_binax_sim: [10 14 20 24 30 34 40 44]
    % hist_sim: [22 24 20 16 8 6 2 2]
    % clustercount_ax: [1×100 double]
    % clustercount_in: [8×1 double]
    % clustercount_out: [8×1 double]
    % clustercount_spread: [8×1 double]
    % clustercount_SEM
 
%examples of input paths
inpth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Paper2018\On cluster_testing\Fig_ClusterTesting Rebuttal 2\V2_FullSetMatlaboutput\';
inpth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Paper2018\On cluster_testing\Fig_ClusterTesting Rebuttal 2\V3_FullSetMatlaboutput\';
inpth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Paper2018\On cluster_testing\Fig_ClusterTesting Rebuttal 2\V4_FullSetMatlaboutput\';
inpth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Matlabcode\2018_09_TestingClusters\Perfectrun\Results_clus_noisefree25reps\';
inpth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Matlabcode\2018_09_TestingClusters\Perfectrun\Results_clus_noiserange25reps\';
inpth='D:\jkerssemakers\_Data_out\2016_Sandro\2018_SJ\2018_09_TestingClusters\Perfectrun\Results_clus_noiserange25reps\';

if 0
%% collect real files
switch 2
    case 1
        realdatapth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Processed\2016_FW\BatchRunData1\BatchanalysisResults\';
        expname='VersionTest';
        realdataexp=['Results_' expname '_Rfp_2Branch_Flipped\'];
    case 2
        realdatapth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Processed\2016_FW\BatchRunData2\BatchanalysisResults\';
        expname='BN2179main';
        realdataexp=['Results_' expname '_Rfp_SingleBranch_NoFlip\']; 
end

realdatafile1=[expname '_A015_Cell_ScreenReport.xlsx']; 
realdatafile2=[expname '_A050_Cell_ClusterReport.xlsx'];

[numdat1,textdat1]= xlsread(strcat(realdatapth,realdataexp,realdatafile1));
[numdat2,textdat2]= xlsread(strcat(realdatapth,realdataexp,realdatafile2));
Headers1=textdat1(1,:);
Headers2=textdat2(2,:);


realclusterno=(numdat2(3:end,(find(strcmp(Headers2,'cluster number')))));
realclusterlabel=(numdat2(3:end,(find(strcmp(Headers2,'label')))));

orilabel=(numdat1(3:end,(find(strcmp(Headers1,'label')))));
realdepth_all=(numdat1(3:end,(find(strcmp(Headers1,'Ratio')))));

LC=length(realclusterlabel);
[~,idx]=ismember(realclusterlabel,orilabel);

realdonutdepth=realdepth_all(idx);
real_avdepth=mean(realdonutdepth);

cnt=0;
for ii=1:10
   sel=find(realclusterno==ii); Lcl=length(sel);
   if Lcl>0
        cnt=cnt+1;
        realclus_no(cnt)=ii;
        realclusdepth_av(cnt)=nanmean(realdonutdepth(sel));
        realclusdepth_st(cnt)=nanstd(realdonutdepth(sel));
        realclusdepth_sem(cnt)=realclusdepth_st(cnt)/(Lcl)^0.5;
   end
end
end

%Simulations  
close all;
switch 1
    case 0
        noizes=[0]
        measnoizes=[0];
        noizlabel=[{'0'}];
        plotax=[1]';
    case 1
        noizes=[0  0.02  0.05  0.2 0.5]
        measnoizes=[0 0.003 0.01 0.04 0.09];
        noizlabel=[{'0'},{'0.003'},{'0.01'},{'0.04'},{'0.09'}];
        plotax=[1 2 3 4 5]';
    case 2
        noizes=[0   0.05  0.5]
        measnoizes=[0  0.01  0.09];
        noizlabel=[{'0'},{'0.01'},{'0.09'}];
        plotax=[1 2 3]';
end
LN=length(noizes);

    


if 0
 %% 1)collect and process data per cell    
    data=load(strcat(inpth,'\all_data_per_cell.mat'),'all_data_per_cell');
    Nc_in=data.all_data_per_cell(:,1);
    Nc_out=data.all_data_per_cell(:,2);
    noiz=data.all_data_per_cell(:,3);
    lomidhi=[0 0.1 0.3 100];
    LR=length(lomidhi)-1;
    clusax=(1:10)';
    errorcurves=NaN*zeros(10,LR);
    for nz=1:LR  %for all noise ranges
        lo=lomidhi(nz); hi=lomidhi(nz+1);
        for ii=1:8  %for all cluster sizes       
            sel=find((ii==Nc_in)&(noiz>=lo)&(noiz<hi));  %noise range per cluster size
            if~isempty(sel)
                errorcurves(ii,nz)=nanmean(abs(Nc_out(sel)));
            end
        end
    end
    subplot(2,2,1);
          plot(clusax,errorcurves,'o-'); hold on;
          legend('error');
          xlabel('cluster number');
          ylabel('error');
          legend('low', 'mid', 'high');
          %axis([0 0.5 0 20]);
          axis square
end


%% 2) collect per input noise
for noi=1:LN;
    inname=strcat('distrib_in_out_noise',num2str(noizes(noi)));
    target_plotdata=strcat(inpth,inname,'_plotdata.mat');
    load(target_plotdata,'Plotdata', 'sim');
    if noi==1
        clustercount_in=repmat(Plotdata.clustercount_in',1,LN);
        [Nc,~]=size(clustercount_in);
        clustercount_out=zeros(Nc,LN);       
        clustercount_sem=zeros(Nc,LN);
        clustercount_spread=zeros(Nc,LN);
        holedepth_av=zeros(Nc,LN);
        holedepth_sem=zeros(Nc,LN);
    end
    clustercount_out(:,noi)=Plotdata.clustercount_out;    
    clustercount_sem(:,noi)=Plotdata.clustercount_SEM;
    clustercount_spread(:,noi)=Plotdata.clustercount_spread;
    holedepth_av(:,noi)=Plotdata.donutness_av;
    holedepth_sem(:,noi)=Plotdata.donutness_SEM;
end
clustercount_dif=clustercount_out-clustercount_in;
failrange=find(abs(clustercount_dif)>0.75);
goodrange=find(abs(clustercount_dif)<=0.75);

holedepth_av_goodrange=holedepth_av;
holedepth_av_goodrange(failrange)=NaN;
holedepth_av_failrange=holedepth_av;
holedepth_av_failrange(goodrange)=NaN;

if 0
binax=-0.02:0.04:0.5;
noisehist=hist(realdonutdepth,binax);
noisehist=noisehist/sum(noisehist)*100;
end

error_for_all_clusters=mean(clustercount_dif)




%% plotting

figure;
 
subplot(2,2,1);
    plotax=[1]'; 
   errorbar(clustercount_in(:,plotax),clustercount_out(:,plotax),clustercount_spread(:,plotax),...
            'o', 'Linewidth', 1, 'Color', [0.5 0.5 0.5]); hold on;
    errorbar(clustercount_in(:,plotax),clustercount_out(:,plotax),clustercount_sem(:,plotax),...
            'ro', 'Linewidth', 1, 'MarkerSize',3); hold on;
    xlabel('Cluster Count In');
    ylabel('Cluster Count Out');
    plot(0:7,(0:7),'k-'); hold on;
    axis square
    axis([1 7 1 7]);

    subplot(4,4,4);
 
    plotax=[2  3  5]';
    errorbar(clustercount_in(:,plotax),clustercount_dif(:,plotax),clustercount_sem(:,plotax),...
            'o-', 'Linewidth', 1, 'MarkerSize',3); hold on;
    xlabel('Cluster Count In');
    ylabel('Counting Error');
    legend(noizlabel(plotax));
    plot(0:7,0*(0:7),'k-'); hold on;
    axis square
    axis([1 7 -1 1]);


subplot(2,2,3);  
    plotax=[1 2 3 4 5]';
    plot(measnoizes,error_for_all_clusters, 'bo-'); hold on;
    plot([0.01 0.01]', [-0.25 0.0]','r--');
    legend('simulations', 'real cells');
    xlabel('Measured Noise');
    ylabel('Average Error');
    ylim([-0.25 0.0]);
    axis square

      
      % subplot(2,2,3);
