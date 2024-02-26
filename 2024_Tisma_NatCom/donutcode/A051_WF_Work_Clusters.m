function A051_WF_Work_Clusters(batchrunindex)
%JWJK_B:-------------------------------------------------------------------
%Plotting of 1D density curves
%
%Summary: This analysis produces heat maps of the collected aligned 
%chromosome maps built by 'A020'and plots the result. further, the curves
%are grouped in 2D 'kymographs' (position-index heat maps). 
%
%Approach: The kymographs have the option to be sorted for various features. 
%Furthermore, curves can be periodically padded to highlight features near
%ori (set in the 'initval' file).
%
%Input: pre-stored .mat data.
%Output: summary plots and excel tables containing the curves.
%
%:JWJK_B-------------------------------------------------------------------
actions.save_svg_wrapup=0; 
actions.save_jpg=0;


initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);
A051_resultpath=[initval.resultpath,initval.DirSep, 'A051_cluster_relations\'];
mkdir(A051_resultpath);   

%% run all combinations:
%2=DNA 3=SMC 4=parB
%'412' means c4 main ori - clusters c2  etc.
combinations=[212 213 223 312 322 402 412 413 423];  %check below 

%% rebuttal work I:
%1st DNA vs 1st SMC: 312
%2nd DNA vs 2nd SMC: 322
combinations=[312 322];  %check below 
%combinations=[412];  %check below
%%
%% panel work

for cii=1:length(combinations)
    thiscombi=combinations(cii);
    %choose how to set up cluster coordination:     
    switch thiscombi
       
        case 312 
            titl='DNA clusters vs SMC-main cluster';
            actions.reference_modus='c3_maincluster';
            actions.clusterplot_modus='c2';        
        case 322 
            titl='DNA clusters vs SMC-secondary cluster';
            actions.reference_modus='c3_secondcluster';
            actions.clusterplot_modus='c2';
        case 212 
            titl='DNA clusters vs DNA-primary cluster';
            actions.reference_modus='c2_maincluster';
            actions.clusterplot_modus='c2';
        case 213 
           titl='SMC clusters vs DNA-primary cluster';
            actions.reference_modus='c2_maincluster';
            actions.clusterplot_modus='c3';
        case 223 
            titl='SMC clusters vs DNA-secondary cluster';
            actions.reference_modus='c2_secondcluster';
            actions.clusterplot_modus='c3';
        case 402 
            titl='DNA clusters vs parB-main spot';
            actions.reference_modus='c4_singlespot'; %use the single-spot parB xy position
            actions.clusterplot_modus='c2';    
        case 412 
            titl='DNA clusters vs parB-main cluster';
            actions.reference_modus='c4_maincluster'; %choose how to set up cluster coordination: 
            actions.clusterplot_modus='c2';       
        case 413 
            titl='SMC clusters vs ParB-primary cluster';
            actions.reference_modus='c4_maincluster';
            actions.clusterplot_modus='c3';
        case 423 
            titl='SMC clusters vs ParB-secondary cluster';
            actions.reference_modus='c4_secondcluster';
            actions.clusterplot_modus='c3';
    end


    %in microns around parB(or equivalent main peak (supposedly centered)
    %We estimate ParB-associated DNA peak at +0.2 micron
    pk=0.2;
    binslots=4;
    binwidth=0.4;
    croplim=[pk-binwidth pk+binslots*binwidth];


    svg_exports=initval.svg_exports;

    close all; pause(0.1);
    initval.Cell_Labels;
    allframes=length(initval.Cell_Labels);
    goodcount=0;
    allmaxima=[];
    cluster_wrapup_bright_nonbright=[];
    cluster_fraction_nearby_far=[];
    contourlengths=[];
    %loop 2: paste spaghetti curves by sorting index
    Nframes=min([200E6 allframes]);
    for jj=1:Nframes
        cellno=char(initval.Cell_Labels{jj});   
        do_this_cell=1;
        if do_this_cell
            disp(strcat('A051_Analyzing cluster relations..cell',num2str(jj),'and',num2str(Nframes-jj+1),'cells to go:'));
            CellName=strcat('ResultsOfCell',cellno,'.mat'); 
            MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
            load(strcat(MatFilePath,CellName));
            if (GeneralCellProps.Okayproduct|initval.PassAllCells)&Aligned.BP.CorrectionOK
                goodcount=goodcount+1;    
                close all;
                switch actions.clusterplot_modus
                    case 'c2', Clusters=Clusters_c2;
                    case 'c3', Clusters=Clusters_c3;
                    case 'c4', Clusters=Clusters_c4;
                end
                [ref_x,ref_y, sub_select_this_cell]=get_ref(c4, Clusters_c2, Clusters_c3, Clusters_c4,actions);
                
                   if sub_select_this_cell
                       %2) construct a genome contour to use as linear coordanate
                        %geometrical center:
                        xCom=Chromosome.xCOM;
                        yCom=Chromosome.yCOM;
                        %relevant contours: 
                        chrx=Chromosome.CartesianContourMax_X;
                        chry=Chromosome.CartesianContourMax_Y; 
                        chr_edge_x=Chromosome.CartesianContourEdge_X;
                        chr_edge_y=Chromosome.CartesianContourEdge_Y;
                        cwx=Chromosome.CartesianCellwallX;
                        cwy=Chromosome.CartesianCellwallY;
                        sel=find(~isnan(cwx));
                        cwx=cwx(sel);
                        cwy=cwy(sel);
                        %measure lengths:    
                        distax0=[0 [cumsum(((chrx(2:end)-chrx(1:end-1)).^2+(chry(2:end)-chry(1:end-1)).^2).^0.5)]'];
                        distax_edge=[0 [cumsum(((chr_edge_x(2:end)-chr_edge_x(1:end-1)).^2+(chr_edge_y(2:end)-chr_edge_y(1:end-1)).^2).^0.5)]'];
                        distax_cw=[0 [cumsum(((cwx(2:end)-cwx(1:end-1)).^2+(cwy(2:end)-cwy(1:end-1)).^2).^0.5)]'];
                        %full lengths of contours:
                        Lc=round(nanmax(distax0));
                        Lc_edge=round(nanmax(distax_edge));
                        Lc_cw=round(nanmax(distax_cw));
                        contourlengths=[contourlengths; ...
                                        [initval.nmperpixel*[Lc Lc_edge Lc_cw]/1000]];

                        %equalize along said contour:
                        [contour_xe,contour_ye]=B002_EqualizeAlongContour(chrx,chry,Lc); 

                        %adapt start position (on non-aligned data):
                        idx=find_contour_start(ref_x,ref_y,contour_xe,contour_ye);  

                        %peak contour (pixel units):
                        px=circshift(contour_xe,-idx);
                        py=circshift(contour_ye,-idx);
                        %distances forward and backward 
                        p_dist_pix_forward=[0 [cumsum(((px(2:end)-px(1:end-1)).^2+(py(2:end)-py(1:end-1)).^2).^0.5)]'];
                        pxb=flipud(px);
                        pyb=flipud(py);
                        p_dist_pix_backward=[0 [cumsum(((pxb(2:end)-pxb(1:end-1)).^2+(pyb(2:end)-pyb(1:end-1)).^2).^0.5)]'];
                        p_dist_pix_backward=fliplr(p_dist_pix_backward);

                        %for plotting purposes: cell wall
                        BX=circshift(Chromosome.CartesianCellwallX,-idx);
                        BY=circshift(Chromosome.CartesianCellwallY,-idx);

                       %to orient the clusters, first map up the cell via the DNA
                       %pattern:

                       %sample dna images:
                       [ICX,ICY,OCX,OCY]=build_extra_contour_rings(px,py,xCom,yCom);
                       c2_map=build_gridmap(c2_pic,ICX,ICY,OCX,OCY);

                       %set reference in center:
                       centershift=round(length(ICX)/2);
                       c2_map=circshift(c2_map,-centershift,2);

                       %assume dark space:
                       c2_map=c2_map-min(c2_map(:));

                       %symmetry check:
                       flipmodus='localsymmetry' ; %'classic';
                       switch flipmodus
                           case 'classic', needtoflip=strcmp(Aligned.Orientation, 'Tails')
                           case 'localsymmetry', needtoflip=check_local_symmetry(c2_map,0.5); 
                       end         
                       if needtoflip  
                            px=flipud(px);
                            py=flipud(py);
                            BX=flipud(BX);
                            BY=flipud(BY);
                       end

                       %link clusters to contour 
                       [~,clust_N]=size(Clusters);
                       clust_x_=zeros(clust_N,1);
                       clust_y=zeros(clust_N,1);
                       clust_perc=zeros(clust_N,1);
                       clust_dist_microns=zeros(clust_N,1);
                       clust_brightest=zeros(clust_N,1);
                       %find nearest backbone position 
                       for ci=1:clust_N
                           clust_x_(ci)=Clusters(ci).COM_X;
                           clust_y(ci)=Clusters(ci).COM_Y;
                           clust_brightest(ci)=ci;
                           clust_perc(ci)=round(100*Clusters(ci).C_perc);
                           idx=find_contour_start( clust_x_(ci), clust_y(ci),px,py);
                           dist_bw=initval.nmperpixel*p_dist_pix_backward(idx)/1000;
                           dist_fw=initval.nmperpixel*p_dist_pix_forward(idx)/1000;
                           %if cluster nearby is on other side, assume periodicity
                           %issue
                           if dist_bw<dist_fw && dist_bw<-croplim(1)
                               clust_dist_microns(ci)=-dist_bw;
                           else
                               clust_dist_microns(ci)=dist_fw;
                           end  
                       end
                       clust_perc=clust_perc/sum(clust_perc)*100;
                       plotscat=0.1*(rand(1,1)-0.5);
                       clust_dist_microns=clust_dist_microns+plotscat;
                       %1 collect Clusters_c2 and indicate brightest
    %                    [~, bright_i]=max(clust_perc);
    %                    clust_brightest(bright_i)=1;
                       cluster_wrapup_bright_nonbright=[cluster_wrapup_bright_nonbright; [clust_dist_microns clust_perc clust_brightest]];

                       %2 collect intensities in nearby and far away, per cell!
                       sel_nearby=find(abs(clust_dist_microns)<=0.5);
                       sel_far=find(abs(clust_dist_microns)>0.5);
                       nearby_fraction=sum(clust_perc(sel_nearby));
                       faraway_fraction=sum(clust_perc(sel_far));
                       sumfraction=nearby_fraction+faraway_fraction;
                       newentry=100*[nearby_fraction faraway_fraction]/sumfraction;
                       cluster_fraction_nearby_far=[cluster_fraction_nearby_far;...
                           newentry];   
                   end
            end
        end
    end


    %% panel: cluster_wrapup:
    %cluster_wrapup=[dist_clust_microns perc_clust brightest]];
    brightest_clusters=cluster_wrapup_bright_nonbright(cluster_wrapup_bright_nonbright(:,3)==1,:);
    secondary_clusters=cluster_wrapup_bright_nonbright(cluster_wrapup_bright_nonbright(:,3)==2,:);
    other_clusters=cluster_wrapup_bright_nonbright(cluster_wrapup_bright_nonbright(:,3)>2,:);
    figure('Units','normalized','Position',[0 0 1 1]);
    subplot(1,2,1);
        plot(other_clusters(:,1),other_clusters(:,2), 'bo','MarkerFaceColor', 'b','MarkerSize', 2); hold on;   
        plot(secondary_clusters(:,1),secondary_clusters(:,2), 'go', 'MarkerFaceColor', 'g','MarkerSize', 3); hold on;
        plot(brightest_clusters(:,1),brightest_clusters(:,2), 'ro', 'MarkerFaceColor', 'r','MarkerSize', 3); hold on;
        title(titl);
        ylabel('cluster content, %'); 
        xlabel('distance, microns');
        legend('weak clusters', 'second brightest clusters', 'brightest clusters', 'Location', 'NorthOutside');

    subplot(2,2,2);
        cluster_fraction_nearby=cluster_fraction_nearby_far(:,1);
        cluster_fraction_far=cluster_fraction_nearby_far(:,2);
        nearby_perc_av=mean(cluster_fraction_nearby);
        nearby_perc_std=std(cluster_fraction_nearby);
        far_perc_av=mean(cluster_fraction_far);
        far_perc_std=std(cluster_fraction_far);
        bar_ax=[1 2];
        bar(bar_ax,[nearby_perc_av far_perc_av]); hold on;
        errorbar([1 2],[nearby_perc_av far_perc_av],[nearby_perc_std far_perc_std]);
        text(0.5, 80,['N=', num2str(length(cluster_fraction_nearby_far(:,1)))]);
        ylabel('pattern compaction'); 
        xlabel('<0.5/>0.5 mu distance of ref');
        %for tabular saving:
        plotpanel_distance_headers=[...
            {'category'},...
            {'bin index'},...
            {'perc'}...
            {'std perc'}...
            ];
        plotpanel_distance_rownames=[{'within 0.5 mu'},{'beyond 0.5 mu'}]'; 
        plotpanel_distance_data=[bar_ax' [nearby_perc_av far_perc_av]' [nearby_perc_std far_perc_std]'];
        

    subplot(2,2,4);
        N_totaal=length(cluster_fraction_nearby);
        N_over40=length(find(cluster_fraction_nearby>40));
        Perc_over_40=N_over40/N_totaal*100;
        Perc_under_40=100-Perc_over_40;
        bar_ax=[1 2];
        bar(bar_ax,[Perc_over_40 Perc_under_40]); hold on;
        text(0.5, 40,['N=', num2str(length(cluster_fraction_nearby_far(:,1)))]);
        text(0.5, 30,['%over40%=', num2str(Perc_over_40)]);
        text(0.5, 20,['%under40%=', num2str(Perc_under_40)]);
        ylabel('cell percentage'); 
        xlabel('>40%/<40% pattern content near refererence');
        %for tabular saving:
        plotpanel_dense_nondense_cellcounts_headers=[...
            {'category'},...
            {'bin index'},...
            {'perc'}...
            ];
        plotpanel_dense_nondense_cellcounts_rownames=[{'> 40% ontent within 0.5 mu'},{'< 40% content within 0.5 mu'}]'; 
        plotpanel_dense_nondense_cellcounts_data=[bar_ax' [Perc_over_40 Perc_under_40]'];
   
    savename=strcat(A051_resultpath,initval.expi,'_A051_',actions.clusterplot_modus,'_clusters_vs_',actions.reference_modus);
    subsavename=[initval.expi,'_A051_',actions.clusterplot_modus,'_clusters_vs_',actions.reference_modus]
    saveas(gcf,[savename,'.jpg'],'jpg');
    if actions.save_svg_wrapup
        plot2svg([savename,'.svg'], gcf);
    end  

    % save bar data as excel
    %per distance, percentage per cell:
    xlswrite([savename,'.xlsx'], plotpanel_distance_headers,'PerDistance','A1');
    xlswrite([savename,'.xlsx'], plotpanel_distance_rownames,'PerDistance','A2');
    xlswrite([savename,'.xlsx'],plotpanel_distance_data,'PerDistance','B2');
    
    %per cell, dense or extended:
    xlswrite([savename,'.xlsx'], plotpanel_dense_nondense_cellcounts_headers,'PerCell','A1');
    xlswrite([savename,'.xlsx'], plotpanel_dense_nondense_cellcounts_rownames,'PerCell','A2');
    xlswrite([savename,'.xlsx'],plotpanel_dense_nondense_cellcounts_data,'PerCell','B2');
    
    %save scatter data:
    panel_brightest_clusters=brightest_clusters(:,1:2);
    panel_secondary_clusters=[secondary_clusters(:,1:2)];
    panel_other_clusters=[other_clusters(:,1:2)];
    
    %per cell, dense or extended:
    scatter_headers=[{'distance'},{'fraction %'}];
    xlswrite([savename,'.xlsx'], scatter_headers,'Scatter_brightests','A1');
    xlswrite([savename,'.xlsx'], panel_brightest_clusters,'Scatter_brightests','A2');
    xlswrite([savename,'.xlsx'], scatter_headers,'Scatter_secondaries','A1');  
    xlswrite([savename,'.xlsx'],panel_secondary_clusters,'Scatter_secondaries','A2');
    xlswrite([savename,'.xlsx'], scatter_headers,'Scatter_others','A1');  
    xlswrite([savename,'.xlsx'],panel_other_clusters,'Scatter_others','A2');
    
     %rebuttal work:
     if 1 & (thiscombi== 312 | thiscombi==322)
         
        %run separate code generation lines or save svg
        close all;
        figure(337);
        subplot(1,3,1);
        plot(secondary_clusters(:,1),secondary_clusters(:,2), 'o', 'MarkerEdgeColor', [0.5 0.5 0.5],'MarkerSize', 2); hold on;
        axis square
 
        ylabel('cluster content, %'); 
        xlabel('distance, microns');
        ylim([0 100]);
        xlim([-0.5 3.5]);
        subsavename=[initval.expi,'_A051_',actions.clusterplot_modus,'_clusters_vs_',actions.reference_modus]
        saveas(gcf,[subsavename,'.jpg'],'jpg');
        %add_a_break_here; 
        dum=1 
     end
    

end

function [ref_x,ref_y, sub_select_this_cell]=get_ref(c4, Clusters_c2, Clusters_c3, Clusters_c4,actions);
%choose origin_of_interest related to a spot postion               
    switch actions.reference_modus % 
        case 'c4_singlespot' %use al approved cells, parB spot leads
            ref_x=c4.spotX; 
            ref_y=c4.spotY;
            sub_select_this_cell=1;
        case 'c4_maincluster' %use al approved cells, parB main cluster leads
            if isfield(Clusters_c4,'COM_X')
                ref_x=Clusters_c4(1).COM_X;
                ref_y=Clusters_c4(1).COM_Y;
                sub_select_this_cell=1;
            else
                ref_x=NaN;
                ref_y=NaN;
                sub_select_this_cell=0;
            end
        case 'c4_secondcluster' %use al approved cells, SMC main cluster leads
            if isfield(Clusters_c4,'COM_X')
                [~,Nc]=size(Clusters_c4);
                if Nc>1
                    ref_x=Clusters_c4(2).COM_X;
                    ref_y=Clusters_c4(2).COM_Y;
                    sub_select_this_cell=1;
                else
                    ref_x=NaN;
                    ref_y=NaN;
                    sub_select_this_cell=0;
                end
            else
                ref_x=NaN;
                ref_y=NaN;
                sub_select_this_cell=0;
            end
        case 'c3_maincluster' %use al approved cells, SMC main cluster leads
            if isfield(Clusters_c3,'COM_X')
                ref_x=Clusters_c3(1).COM_X;
                ref_y=Clusters_c3(1).COM_Y;
                sub_select_this_cell=1;
            else
                ref_x=NaN;
                ref_y=NaN;
                sub_select_this_cell=0;
            end
        case 'c3_secondcluster' %use al approved cells, SMC main cluster leads
            if isfield(Clusters_c3,'COM_X')
                [~,Nc]=size(Clusters_c3);
                if Nc>1
                    ref_x=Clusters_c3(2).COM_X;
                    ref_y=Clusters_c3(2).COM_Y;
                    sub_select_this_cell=1;
                else
                    ref_x=NaN;
                    ref_y=NaN;
                    sub_select_this_cell=0;
                end
            else
                ref_x=NaN;
                ref_y=NaN;
                sub_select_this_cell=0;
            end
        case 'c2_maincluster' %use al approved cells, parB main cluster leads
            if isfield(Clusters_c2,'COM_X')
                ref_x=Clusters_c2(1).COM_X;
                ref_y=Clusters_c2(1).COM_Y;
                sub_select_this_cell=1;
            else
                ref_x=NaN;
                ref_y=NaN;
                sub_select_this_cell=0;
            end
        case 'c2_secondcluster' %use al approved cells, SMC main cluster leads
            if isfield(Clusters_c2,'COM_X')
                [~,Nc]=size(Clusters_c2);
                if Nc>1
                    ref_x=Clusters_c2(2).COM_X;
                    ref_y=Clusters_c2(2).COM_Y;
                    sub_select_this_cell=1;
                else
                    ref_x=NaN;
                    ref_y=NaN;
                    sub_select_this_cell=0;
                end
            else
                ref_x=NaN;
                ref_y=NaN;
                sub_select_this_cell=0;
            end
    end 
     

           