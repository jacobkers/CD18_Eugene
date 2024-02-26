function A055_WF_Work_Demographs(batchrunindex)
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

actions.workonlyexamplecells=0;
actions.save_svg_wrapup=1; %only if init is on
actions.save_jpg=1;
actions.show_also_SMC=1;

if ~actions.workonlyexamplecells
    selectlabel='_all';
else
    selectlabel='_selected';
end


 %in microns around parB (supposedly centered)
%We estimate ParB-associated DNA peak at +0.2 micron
pk=0.2;
binslots=4;
binwidth=0.4;
croplim=[pk-binwidth pk+binslots*binwidth];


selectcurves_idx=[];
selectcurves_c2=[];
greymap=repmat([0:0.1:1.0]',1,3);
greenmap=greymap;
greenmap(:,1)=0;
greenmap(:,3)=0;

orangemap=greymap;
orangemap(:,1)=239/255*orangemap(:,1);
orangemap(:,2)=156/255*orangemap(:,2);
orangemap(:,3)=0;


initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);

%clustersource=(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A050_ClusterGeneralprops.mat'));


imoutdir=strcat(initval.resultpath,'A35_CellImages_Crescent_Analysis',initval.DirSep);   
%if isdir(imoutdir), rmdir(imoutdir,'s');  end
exampledir=[imoutdir, '\ExampleCells\'];
if ~isdir(imoutdir), mkdir(imoutdir);  end
if ~isdir(exampledir), mkdir(exampledir);  end
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
    %check if it is an examplecell:
    use_as_example=~isempty(find(strcmp(cellno,svg_exports)==1));    
    do_this_cell=((actions.workonlyexamplecells&use_as_example)|...
                 ~actions.workonlyexamplecells);
    if do_this_cell
        disp(strcat('Building demographs:cell',num2str(jj),'and',num2str(Nframes-jj+1),'cells to go:'));

        CellName=strcat('ResultsOfCell',cellno,'.mat'); 
        MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
        load(strcat(MatFilePath,CellName));               
       %% Initialize and fill two collections of chromosome density curves:   

        if (GeneralCellProps.Okayproduct|initval.PassAllCells)&Aligned.BP.CorrectionOK
                goodcount=goodcount+1;        
                if goodcount==1     
                    laa=100;
                    demograph_c2=zeros(Nframes,laa);
                    demograph_c3=zeros(Nframes,laa);
                    demograph_c4=zeros(Nframes,laa);
                end   

               close all;
               switch initval.alignmodus_A55
                   case 'c4'
                        spot_alignx=c4.spotX; 
                        spot_aligny=c4.spotY;
                   case 'c3'
                        spot_alignx=c3.spotX; 
                        spot_aligny=c3.spotY;
               end
               
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
                idx=find_contour_start(spot_alignx,spot_aligny,contour_xe,contour_ye);  
               
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
               
               
               %% 1D-mapping section 
               %set full-round contours:
               [ICX,ICY,OCX,OCY]=build_extra_contour_rings(px,py,xCom,yCom);

               %sample images:
               phase_map=build_gridmap(cell_pic,ICX,ICY,OCX,OCY);
               c2_map=build_gridmap(c2_pic,ICX,ICY,OCX,OCY);
               c3_map=build_gridmap(c3_pic,ICX,ICY,OCX,OCY);
               c4_map=build_gridmap(c4_pic,ICX,ICY,OCX,OCY);
                              
               %sample images:
               c2_map=build_gridmap(c2_pic,ICX,ICY,OCX,OCY);
               c3_map=build_gridmap(c3_pic,ICX,ICY,OCX,OCY);
               c4_map=build_gridmap(c4_pic,ICX,ICY,OCX,OCY);

               %set parB in center:
               centershift=round(length(ICX)/2);
               c2_map=circshift(c2_map,-centershift,2);
               c4_map=circshift(c4_map,-centershift,2);
               c3_map=circshift(c3_map,-centershift,2);
               phase_map=circshift(phase_map,-centershift,2);

               %assume dark space:
               c2_map=c2_map-min(c2_map(:));
               c4_map=c4_map-min(c4_map(:));
               c3_map=c3_map-min(c3_map(:));
               phase_map=phase_map-min(phase_map(:));


               %adapt direction:
               flipmodus='localsymmetry' ; %'classic';
               switch flipmodus
                   case 'classic', needtoflip=strcmp(Aligned.Orientation, 'Tails')
                   case 'localsymmetry', needtoflip=check_local_symmetry(c2_map,0.5); 
               end         

               if needtoflip  %flip order of contours!
                    fliptext='flipped';
                    px=flipud(px);
                    py=flipud(py);
                    BX=flipud(BX);
                    BY=flipud(BY);
                    c2_map=fliplr(c2_map);
                    c4_map=fliplr(c4_map);
                    c3_map=fliplr(c3_map);
                    phase_map=fliplr(phase_map);
               else
                   fliptext='';
               end

               % finalization translate to microns.
               %2) define a distance axis (note we have 'floating' units):
               pixelsperunit=Chromosome.TotalMaxPeakLength/length(px);
               scalefactor=initval.nmperpixel*pixelsperunit/1000;

               %Note: we want to avoid periodic profile contributions
               pts=100;
               [phase_map, ~]=resize_map(phase_map, croplim, scalefactor,pts);
               [c2_map, distax]=resize_map(c2_map, croplim, scalefactor, pts);
               [c4_map, ~]=resize_map(c4_map, croplim, scalefactor,pts);               
               [c3_map, ~]=resize_map(c3_map, croplim, scalefactor,pts);

               %finally, obtain the profiles:
               %DNA:
               profile_to_use_c2=max(c2_map, [], 'omitnan');
               profile_to_use_c2=medsmooth(profile_to_use_c2,2, 'Average');
               profile_to_use_c2=profile_to_use_c2/sum(profile_to_use_c2)*100; 
               %parB
               profile_to_use_c4=max(c4_map, [], 'omitnan');
               profile_to_use_c4=medsmooth(profile_to_use_c4,2, 'Average');
               profile_to_use_c4=profile_to_use_c4/sum(profile_to_use_c4)*100; 
               %SMC:
               profile_to_use_c3=max(c3_map, [], 'omitnan');
               profile_to_use_c3=medsmooth(profile_to_use_c3,2, 'Average');
               profile_to_use_c3=profile_to_use_c3/sum(profile_to_use_c3)*100; 
                           
               maxima=Eval_Maxima(distax,profile_to_use_c2,distax);
               allmaxima=[allmaxima;  [maxima.Idx maxima.Pos maxima.Vals]];
               
               demograph_c2(goodcount,:)=profile_to_use_c2;
               demograph_c3(goodcount,:)=profile_to_use_c3;
               demograph_c4(goodcount,:)=profile_to_use_c4;
                             
               
               %% Plot-per_cell
               %1D-mapping plot menu: 
               if actions.save_jpg 
                %set size
                %switch datamodus
                   [ori_rws,~]=size(c2_map);
                   combi_map=[1000*c4_map/max(c4_map(:))  
                              1000*c2_map/max(c2_map(:))];
                   
                   if actions.show_also_SMC
                       combi_map=[1000*c3_map/max(c3_map(:))  
                                 combi_map]; 
                   end
     
                   [rm,cm]=size(combi_map);
                   combi_N=round(rm/ori_rws);
                   [Xm,Ym]=meshgrid(1:cm,1:rm);
                   
                   
                   
                  figure('Units','normalized','Position',[0 0 1 0.3], 'Visible', 'off');
                  %figure('Units','normalized','Position',[0 0 1 0.3]);    

                    subplot(1,2,2);
                        scalef=(croplim(2)-croplim(1))/pts;
                        midc=(0-croplim(1))/(croplim(2)-croplim(1));
                        Xms=translate_ax(Xm,scalef,midc*cm);
                        Yms=Ym;
                        pcolor(Xms,Yms,combi_map); shading flat; axis tight; hold on;               
                        offset_chro=rm-(combi_N-2)*ori_rws;
                        offset_parB=rm-(combi_N-1)*ori_rws;
                        offset_SMC= rm-(combi_N)*ori_rws;
                    
                        plot_ampli=0.9*(max(Yms(:,1))-min(Yms(:,1)))/combi_N;
                        profile_to_plot_chro=plot_ampli*profile_to_use_c2/max(profile_to_use_c2)+offset_chro;
                        profile_to_plot_parB=plot_ampli*profile_to_use_c4/max(profile_to_use_c4)+offset_parB;
                        profile_to_plot_SMC=plot_ampli*profile_to_use_c3/max(profile_to_use_c3)+offset_SMC;
                     
                        plot(distax,profile_to_plot_chro, 'w-', 'LineWidth',1); hold on;
                        plot(distax,profile_to_plot_parB, 'w-', 'LineWidth',1); hold on;
                        if  actions.show_also_SMC
                            plot(distax,profile_to_plot_SMC, 'w-', 'LineWidth',1); hold on;
                        end
                    %text(distax(1)+0.2, offset_parB+1,'parB', 'Color','white');
                    xlabel('distance from parB focus, microns');
                    %ylabel('sampling units, a.u.');

                    saveas(gcf,strcat(imoutdir,'A35_crescent_map_aligned',initval.alignmodus_A55,'cell',num2str(cellno,'% 3.0f'),'.jpg')); 

                    %save material for example cells:
                    if 0 %use_as_example==1
                        selectcurves_idx=[selectcurves_idx; jj];
                        selectcurves_c2=[selectcurves_c2; profile_to_use_c2'];
                        close(gcf); 
                        figure('Units','normalized','Position',[0 0 1 0.3]);
                        subplot(1,2,1)
                            pcolor(plot_pic); shading flat; axis off; axis equal; axis tight; hold on; 
                            plot(px,py,'y-','LineWidth',1);
                            plot(BX,BY,'w-', 'LineWidth',1); hold on;
                            plot(spot_alignx,spot_aligny, 'ro','MarkerSize', 6, 'MarkerFaceColor', 'r');
                            legend('DNA','backbone','ParB', 'Location','Eastoutside');
                            hold off;   
                            colormap(greenmap);
                            saveas(gcf,strcat(exampledir,'A35_crescent_map__aligned',initval.alignmodus_A55,'cell', num2str(cellno,'% 3.0f'),'_panel_1.jpg')); 
                            if initval.save_svg
                                plot2svg(strcat(exampledir,'A35_crescent_map__aligned',initval.alignmodus_A55,'cell', num2str(cellno,'% 3.0f'),'_panel_1.svg'), gcf);
                                disp('(exported svg_panel 1...)');
                            end
                            close(gcf); 
                            
                            figure('Units','normalized','Position',[0 0 1 0.3]);
                            subplot(1,2,2)
                            pcolor(Xms,Yms,combi_map); shading flat; axis tight; hold on;
                            xlim(croplim);
                            offset_chro=0;
                            offset_parB=min(Yms(:,1));
                            plot_ampli=0.9*(max(Yms(:,1))-min(Yms(:,1)))/2;
                            profile_to_plot_chro=plot_ampli*profile_to_use_c2/max(profile_to_use_c2)+offset_chro;
                            profile_to_plot_parB=plot_ampli*profile_to_use_c4/max(profile_to_use_c4)+offset_parB;
                            plot(distax,profile_to_plot_chro', 'w-', 'LineWidth',1); hold on;
                            plot(distax,profile_to_plot_parB', 'w-', 'LineWidth',1); hold on;
                            xlabel('distance from parB focus, microns');
                            ylabel('radial sampling units, a.u.');
                            hold off;  
                            pause(0.2);
                            
                            %save in two colors:
                            colormap(greenmap);
                            saveas(gcf,strcat(exampledir,'A35_crescent_map_aligned',initval.alignmodus_A55,'cell', num2str(cellno,'% 3.0f'),'_panel_2_green.jpg')); 
                            if initval.save_svg
                                plot2svg(strcat(exampledir,'A35_crescent_map_aligned',initval.alignmodus_A55,'cell', num2str(cellno,'% 3.0f'),'_panel_2_green.svg'), gcf);
                                disp('(exported svg_panel 2...)');
                            end
                            pause(0.2);
                            colormap(orangemap);
                            saveas(gcf,strcat(exampledir,'A35_crescent_map_aligned',initval.alignmodus_A55,'cell', num2str(cellno,'% 3.0f'),'_panel_2_orange.jpg')); 
                            if initval.save_svg
                                plot2svg(strcat(exampledir,'A35_crescent_map_aligned',initval.alignmodus_A55,'cell', num2str(cellno,'% 3.0f'),'_panel_2_orange.svg'), gcf);
                                disp('(exported svg_panel 2...)');
                            end
                            pause(0.2);
                            close(gcf); 
                    end    
                    close(gcf); 
               end
        end
    end
end

%% remove empty (rejected) cells
okcells=(sum(demograph_c2'))>0;
demograph_c2=demograph_c2(okcells,:);
demograph_c3=demograph_c3(okcells,:);
demograph_c4=demograph_c4(okcells,:);
av_curve_c2=mean(demograph_c2, 'omitnan');
av_curve_c3=mean(demograph_c3, 'omitnan');
av_curve_c4=mean(demograph_c4, 'omitnan');
contourlengths=contourlengths(okcells,:);

[demograph_c2_sort, demograph_c3_sort, demograph_c4_sort]=sort_demographs(demograph_c2, demograph_c3,demograph_c4, 'default');
[rk,ck]=size(demograph_c2_sort);

%% collect random example curves
N_random=20;
if isempty(selectcurves_idx) 
    selectcurves_idx=ceil(rk*rand(N_random,1));
    selectcurves_c2=demograph_c2(selectcurves_idx,:);
    selectcurves_c3=demograph_c3(selectcurves_idx,:);
    selectcurves_c4=demograph_c4(selectcurves_idx,:);
end

%% divide curves in 'slots'
binnedmap_c2=zeros(rk,binslots);
binnedmap_c3=zeros(rk,binslots);
binnedmap_c4=zeros(rk,binslots);
stepsize=round(ck/binslots);
for ii=1:binslots
        ixlo=1+(stepsize*(ii-1));
        ixhi=ixlo+stepsize-1;
    binned_distax(ii)=mean(distax(ixlo:ixhi));
    for jj=1:rk
        prf_c2=100*demograph_c2(jj,:)/sum(demograph_c2(jj,:)); 
        prf_c3=100*demograph_c3(jj,:)/sum(demograph_c3(jj,:));
        prf_c4=100*demograph_c4(jj,:)/sum(demograph_c4(jj,:));
        binnedmap_c2(jj,ii)=sum(prf_c2(ixlo:ixhi));
        binnedmap_c3(jj,ii)=sum(prf_c3(ixlo:ixhi));
        binnedmap_c4(jj,ii)=sum(prf_c4(ixlo:ixhi));
    end
end

slots_sum_av_c2=mean(binnedmap_c2);
slots_std_av_c2=std(binnedmap_c2);
slots_sum_av_c4=mean(binnedmap_c4);
slots_std_av_c4=std(binnedmap_c4);
slots_sum_av_c3=mean(binnedmap_c3);
slots_std_av_c3=std(binnedmap_c3);


%% Saving and plotting results: 

%% plot panel: averages
figure(365);  

subplot(2,3,1);   
    plot(distax,selectcurves_c2,'-','LineWidth',1, 'Color', [0.6 0.6 0.6]); hold on;   
    plot(distax,av_curve_c2,'k','LineWidth',3); hold on;
    title(initval.channelnames{3});
    axis tight;
    %ylim([0.5 1.5]);
subplot(2,3,2);
    plot(distax,selectcurves_c3,'-','LineWidth',1, 'Color', [0.6 0.6 0.6]); hold on;   
    plot(distax,av_curve_c2,'k','LineWidth',3); hold on;
    plot(distax,av_curve_c3,'g-','LineWidth',3); hold on;   
    title(initval.channelnames{4});
    axis tight;  
    %ylim([0.4 2.1]);
subplot(2,3,3);
    plot(distax,selectcurves_c4,'-','LineWidth',1, 'Color', [0.6 0.6 0.6]); hold on;   
    plot(distax,av_curve_c2,'k','LineWidth',3); hold on;
    plot(distax,av_curve_c3,'g-','LineWidth',3); hold on; 
    plot(distax,av_curve_c4,'y-','LineWidth',3); hold on;
    title(initval.channelnames{5});
    axis tight;
    %ylim([0.2 4]);
 
  % demographs 3x 
    skipit=1;
    demograph_c2_srt_crp=demograph_c2_sort(1:skipit:end,:);
    demograph_c3_srt_crp=demograph_c3_sort(1:skipit:end,:);
    demograph_c4_srt_crp=demograph_c4_sort(1:skipit:end,:);
    [r_cr,c_cr]=size(demograph_c2_srt_crp);
    scalef=(croplim(2)-croplim(1))/pts;
    midc=(0-croplim(1))/(croplim(2)-croplim(1));
    [Xk,Yk]=meshgrid(1:c_cr,1:r_cr);
    Xks=translate_ax(Xk,scalef,midc*c_cr);
    Yks=Yk;
    
subplot(2,3,4);
    pcolor(Xks,Yks, demograph_c2_srt_crp); 
    colormap jet; shading flat;hold on;
    colorbar('NorthOutside')
    ylabel('Cell Number'); 
    xlabel('Distance, Microns');
subplot(2,3,5);
    pcolor(Xks,Yks, demograph_c3_srt_crp); 
    colormap jet; shading flat;hold on;
    colorbar('NorthOutside')
    ylabel('Cell Number'); 
    xlabel('Distance, Microns');
subplot(2,3,6);
    pcolor(Xks,Yks, demograph_c4_srt_crp); 
    colormap jet; shading flat;hold on;
    colorbar('NorthOutside')
    ylabel('Cell Number'); 
    xlabel('Distance, Microns');

    
saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A055_Spaghetti_and_Demographs_aligned',initval.alignmodus_A55,selectlabel,'.jpg'),'jpg');
if actions.save_svg_wrapup && initval.save_svg
    plot2svg(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A055_Spaghetti_and_Demographs_aligned',initval.alignmodus_A55,selectlabel,'.svg'), gcf);
end    
    
 
figure(470);    
subplot(3,2,2);    
    bar(binned_distax,slots_sum_av_c2,'w'); hold on;
    errorbar(binned_distax,slots_sum_av_c2,slots_std_av_c2);
    legend('c2');
subplot(3,2,4);       
    bar(binned_distax,slots_sum_av_c4,'w'); hold on;
    errorbar(binned_distax,slots_sum_av_c4,slots_std_av_c4);
    legend('c4');
 subplot(3,2,6);    
    bar(binned_distax,slots_sum_av_c3,'w'); hold on;
    errorbar(binned_distax,slots_sum_av_c3,slots_std_av_c3);
    legend('c3');


%% plot panel: scatter plots
%figure(442);
figure('Units','normalized','Position',[0 0 0.9 0.9]);
subplot(1,3,1); 
    dna_vs_parb_scatterdata=make_scatplot(binnedmap_c4,binnedmap_c2,binnedmap_c3,'ParB','DNA','SMC',30, 'gray');
subplot(1,3,2); 
    smc_vs_dna_scatterdata=make_scatplot(binnedmap_c2,binnedmap_c3,binnedmap_c4,'DNA','SMC','ParB',40,'gray');
%subplot(1,3,3); 
%    make_scatplot(binnedmap_c4,binnedmap_c3,binnedmap_c2,'ParB','SMC','DNA',30,'gray');
%jpg:
saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A055_scatterplots_3x_aligned',initval.alignmodus_A55,selectlabel,'.jpg'),'jpg');
xlsname_scat=strcat(initval.resultpath,initval.DirSep,initval.expi,'_A055_scatterdata',initval.alignmodus_A55,selectlabel,'.xlsx');
%xls, panel A:
hdrs_scat=[{'ParB-percentage'}, {'DNA-percentage'}, {'marker_size'}];
xlswrite(xlsname_scat,hdrs_scat, 'DNA_vs_ParB', 'A1');
xlswrite(xlsname_scat,dna_vs_parb_scatterdata, 'DNA_vs_ParB', 'A2');
%xls, panel B:
hdrs_scat=[{'DNA-percentage'}, {'SMC-percentage'}, {'marker_size'}];
xlswrite(xlsname_scat,hdrs_scat, 'SMC_vs_DNA', 'A1');
xlswrite(xlsname_scat,smc_vs_dna_scatterdata, 'SMC_vs_DNA', 'A2');
%svg:
if actions.save_svg_wrapup & initval.save_svg
    plot2svg(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A055_scatterplots_demographs3x_aligned',initval.alignmodus_A55,selectlabel,'.svg'), gcf);
end


%% panel: contourlenghts
figure(558); 
subplot(1,2,1);
plot(contourlengths(:,3),contourlengths(:,1), 'ro', 'MarkerFaceColor', 'r','MarkerSize', 3); hold on;
plot(contourlengths(:,3),contourlengths(:,2), 'bo', 'MarkerFaceColor', 'r','MarkerSize', 3); 
plot([0 10],[0 10], 'k-'); hold on;
ylabel('length chromosome circle contour, microns'); 
xlabel('length cell wall, microns');
legend('backbone', 'outer edge', '1:1','Location', 'NorthOutside');
saveas(gcf,strcat(initval.resultpath,initval.DirSep,initval.expi,'_A055_contour_lengths_aligned',initval.alignmodus_A55,selectlabel,'.jpg'),'jpg');
if actions.save_svg_wrapup & initval.save_svg
    plot2svg(strcat(initval.resultpath,initval.DirSep,initval.expi,'_A055_contour_lengths_aligned',initval.alignmodus_A55,selectlabel,'.svg'), gcf);

end 
xlsname=strcat(initval.resultpath,initval.DirSep,initval.expi,'_A055_contour_lengths_aligned',initval.alignmodus_A55,selectlabel,'.xlsx');
hdrs=[{'back bone' {'outer edge'} {'cell wall'}}];
xlswrite(xlsname,hdrs, 'Sheet1', 'A1');
xlswrite(xlsname,contourlengths, 'Sheet1', 'A2');
disp('done');
        
 function [demograph_c2_sort, demograph_c3_sort, demograph_c4_sort]=sort_demographs(demograph_c2, demograph_c3,demograph_c4, theme);
 switch theme
    case   'default'
        %% sort the maps on DNA peak value
        [demograph_c2_sort,sort_idx]=B050_ProcessFurther_Sort_CellMaps(demograph_c2, 'Peak');
        demograph_c3_sort=demograph_c3(sort_idx,:);
        demograph_c4_sort=demograph_c4(sort_idx,:);       
     case   'run104:120231123_BSG5522_DAPI_2_sort_on_ter:Fig.S10'
        %% sort the maps on c4 (ter)
        [demograph_c4_sort,sort_idx]=B050_ProcessFurther_Sort_CellMaps(demograph_c4, 'PeakPos');
        demograph_c3_sort=demograph_c3(sort_idx,:);
        demograph_c2_sort=demograph_c2(sort_idx,:);       
 end
 
     

           