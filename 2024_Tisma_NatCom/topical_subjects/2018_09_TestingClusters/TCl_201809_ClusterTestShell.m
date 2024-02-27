function TCl_201809_ClusterTestShell
%This shell generates a series of cluster-donut patterns and analyzes them
%and plots the output cluster number vs. the input cluster number
%Jacob Kers 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% to use either simulations directly, or have these first deconvolved via
% external programs, set these:
actions.generate=1;     %leave at 1
actions.loadnotmake=1;  %1:files will be loaded from a specified directory
                        %2: files will be simulated and saved
actions.analyze=1;      %perform cluster analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%set noise parameter for simulations; when running from loaded files, these
%should match those specified in the file names
noizes=[0  0.02  0.05  0.2 0.5]; LN=length(noizes);

all_data_per_cell=[];  %to use for collecting information
for noi=1:LN %for all noise levels
close all;
%fixed parameters for creating clusters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%this is copied from the 'generate function, but here the 'sim.noise'
%parameter will be overwritten with the above numbers used in the loop

    %simulation settings (take care that with pre-loaded files, the settings match 
    %those in the filenames
    sim.style='randomized';'randomized';%'randomized';  
    sim.donutwidth=0.6;
    sim.holesize=0.7;  %fill donut 0-1
    sim.structurepsf=4;
    sim.repeat=25;
    sim.noise=0.0;
    sim.noisemode='pixelnoise'; %'extrablobs';  % 'pixelnoise'
    sim.noisepeaksno=100;
    sim.c0=80;            %columns picture
    sim.r0=80;            %rows,in pixels
    sim.blowup=7;          %To suppress sampling bias
    sim.psf=2.5;            %point spread function to blur
    sim.fotons=10000;       %average fotons per spot
    sim.countsperfoton=380;  
    sim.s1=4*sim.psf/sim.structurepsf; %cluster spacing, in units of (strucutral) psf.
    sim.n1=5;   %number of clusters.
    sim.s2=0.05; %in-cluster spot spacing, in units of psf
    sim.n2=round(sim.c0/5/sim.n1/sim.s2*3.5/sim.structurepsf);   %number of spots per cluster. 

    %cluster analysis
    initval.Psf_est=2.7; %2.7;
    initval.skipsmallspots=2/100;  %max fraction of spot to reject
    initval.Separation=4.5; %2.7; 
    initval.outpath=pwd;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %cluster loop: here a range of cluster numbers is simulated or loaded and analyzed    
    sim.noise=noizes(noi);    
    repeats=sim.repeat;
    counter=0;
    for nc_i=2:1:6 %2:1:6
        counter=counter+1;
        ClusterGenprops=struct('cell_label',[]);
        SimClusterContents=[];
        MeasClusterContents=[];
        Donutness=[];
        for rep_i=1:repeats
            disp(['nc_i',num2str(nc_i),'Rep',num2str(rep_i)]);
            sim.repeat=rep_i;
            sim.n1=nc_i;   %number of clusters.     
            if actions.generate                
                 if ~actions.loadnotmake
                    [CellName,im,clu_simprops]=TCl_201809_build_cluster_testpattern(sim);
                 else  %make it
                    [CellName,im,clu_simprops]=TCl_201809_load_cluster_testpattern(sim);
                 end
                SimClusterContents=[SimClusterContents ; clu_simprops.perc];
            end
            initval.CellName=CellName;             
            if actions.analyze
                %do cluster analysis%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                [Nout,AllClusters,Clusters,GenProps]=TCl_201809_get_clusters(im,initval);
                relperc=AllClusters(:,1)/sum(AllClusters(:,1))*100;
                MeasClusterContents=[MeasClusterContents;relperc];
                
                %collect cluster-out distribution%%%%%%%%%%%%%%%%%%%%%%%%%
                ClusterGenprops=TCl_201809_get_clustergenprops(Clusters,GenProps,ClusterGenprops,rep_i);
                Donutness(rep_i)=clu_simprops.donutness;
                Nc_out(rep_i)=Nout;
                dev_C(rep_i)=100*mean(abs(AllClusters(:,1)-1/nc_i));  
            end
        end
        if actions.analyze
        Nc_in(counter)=nc_i;
        Nc_out_av(counter)=nanmean(Nc_out); 
        dev_C_av(counter)=nanmean(dev_C); 
        Nc_st(counter)=nanstd(Nc_out); 
        dev_C_st(counter)=nanstd(dev_C);  
        donutness_av(counter)=nanmean(Donutness);
        donutness_st(counter)=nanstd(Donutness);
        
        Nin_axis=zeros(repeats,1)+nc_i;
        all_data_per_cell=[all_data_per_cell;...
                       [Nin_axis Nc_out' Donutness']];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %optional: plot example histograms of set and measured 
        %cluster numbers halfway the loop
        if 0&&(nc_i==4)
            figure(2);            
            binax=0:2:100;
            hist_sim=hist(SimClusterContents,binax);
            hist_sim=hist_sim/sum(hist_sim)*100;
            hist_meas=hist(MeasClusterContents,binax);
            hist_meas=hist_meas/sum(hist_meas)*100;
            if ~actions.loadnotmake
                subplot(2,3,2); 
                stem(binax(hist_sim>0),hist_sim(hist_sim>0),'bo'); hold on;
                bar(binax,hist_meas, 'r'); hold on;
                stem(binax(hist_sim>0),hist_sim(hist_sim>0),'bo');
                xlabel('Cluster Content (%)')
                ylabel('Frequency (%)');
                legend('Input', 'Output');
                ylim([0 1.2*max(hist_sim)]);
                xlim([0 60]);
                axis square; hold on;
            end
            Plotdata.hist_binax=binax;
            Plotdata.hist_meas=hist_meas;
            Plotdata.hist_binax_sim=binax(hist_sim>0);
            Plotdata.hist_sim=hist_sim(hist_sim>0);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     end
    end
    if actions.analyze
        Nc_SEM=Nc_st/(repeats^0.5);
        donutness_SEM=donutness_st/(repeats^0.5);    
        figure(2);           
        subplot(2,3,1);
            axz=1:100;
            plot(axz,axz,'k-'); hold on;
            errorbar(Nc_in,Nc_out_av,Nc_st,...
                    'o', 'Linewidth', 1, 'Color', [0.5 0.5 0.5]); hold on;
            errorbar(Nc_in,Nc_out_av, Nc_SEM,...
                    'ro', 'Linewidth', 1, 'MarkerSize',3); hold on;
            xlabel('Cluster Count In');
            ylabel('Cluster Count Out');
            axis square
            axis([0 9 0 9]);       
        if ~actions.loadnotmake
        subplot(2,3,3);
            axz=1:100;
            plot(axz,axz,'k-'); hold on;
            errorbar(Nc_in,donutness_av,donutness_st,...
                    'o-', 'Linewidth', 1, 'Color', [0.5 0.5 0.5]); hold on;
            errorbar(Nc_in,donutness_av,donutness_SEM,...
                    'r-', 'Linewidth', 1); hold on;
            xlabel('Cluster Count In');
            ylabel('donut depth, a.u.');
            axis square
            axis([1 9 0 1]);    
        end 
        
        %collect data for plotting via 'noiseplotter' program
        Plotdata.clustercount_ax=axz; 
        Plotdata.clustercount_in=Nc_in;  
        Plotdata.clustercount_out=Nc_out_av; 
        Plotdata.clustercount_spread=Nc_st; 
        Plotdata.clustercount_SEM=Nc_SEM; 
        Plotdata.donutness_av=donutness_av;    
        Plotdata.donutness_st=donutness_st; 
        Plotdata.donutness_SEM=donutness_SEM; 
        
        %saving section for this noise level  %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        outname=strcat('distrib_in_out_noise',num2str(sim.noise));
        outpth='D:\jkerssemakers\_Recent\CD\BN_CD16_Sandro\Matlabcode\2018_09_TestingClusters';
        target_pic=strcat(outpth,'\Results_clus\',outname,'.jpg');
        saveas(gcf,target_pic);
        target_plotdata=strcat(outpth,'\Results_clus\',outname,'_plotdata.mat');
        save(target_plotdata,'Plotdata', 'sim');
    end   
end 

%saving for all loops%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if actions.analyze
    target_alldata=strcat(outpth,'\Results_clus\all_data_per_cell.mat');
    save(target_alldata,'all_data_per_cell');
end
