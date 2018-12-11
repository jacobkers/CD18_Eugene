
function A020_KymoStackTrack
%JWJK_A:-------------------------------------------------------------------
%Title:A020_KymoStackTrack
%c
%Summary: %This function builds a kymograph from images 
%and analyzes it
%JacobKers2018
%
%:JWJK_A-------------------------------------------------------------------
actions.analyze=1; 
    %0 load pre-analyzed data for plot menu
    %1 re-analyze
actions.plotmodus=2;  
    %0 don't plot or save, 
    %1 save plot 1, don't show
    %2 save plot 2, don't show

initval.psf_est=3;
    
close all;
loadImageJ_kymograph=1;
datapath='D:\jkerssemakers\_Data\CD\2018_Eugene\';
for exprun=1:3;
switch exprun
    case 1
        expname='2018_08_01 Pilot Runs';       
        AllExp=[1 3];  %numbers of various rois   
    case 2
        expname='2018_09_24 More_molecules';
        AllExp=[2 4 5 6 7 8 9 10 11 12 13 14 15 16];  %paths to various rois  
        %AllExp=[7 8 9 10 11 12 13 14 15 16];  %paths to various rois   
        %AllExp=[14];  %paths to various rois 
    case 3
        expname='2018_12_05 EvenMoreMolecules';
        AllExp=[21 22 23 24 25 26 27];  %paths to various rois  
end
generaldatapth=[datapath,expname,'\'];
outpath=strcat(datapath, 'matlabresults\',expname,'\');
 

if ~isdir(outpath), mkdir(outpath); end

%% standardized subdirectory names 
Channel_list=[{'DNA\'}, {'Condensin\'}];     %The two subchannels
Kymo_list='kymo_ImageJ\'; %ImageJ-made kymographs
Condensin_Kymo='Kymograph_Condensin.tif';           %if you use it
Dna_Kymo='Kymograph_DNA.tif';                %if you use it';


%% main loop
LE=length(AllExp);  %for all experiments
if actions.analyze
for ee=1:LE
if mod(ee,1)==0, disp(strcat('Analyzing:',expname,':Exps to work through:',num2str(LE-ee+1)));end 
Exp=strcat('ROI',num2str(AllExp(ee)));
SaveName=char(strcat(outpath, Exp));  
% 1) build stack of images and kymograph (or load it)
if loadImageJ_kymograph             %load kymograph
    dna_name=char(strcat(generaldatapth, Exp, '\', Kymo_list,Dna_Kymo));
    condensin_name=char(strcat(generaldatapth, Exp,'\', Kymo_list,Condensin_Kymo));
    if exist(condensin_name)==2, do_condensin=1;, else  do_condensin=0;end   
    if do_condensin, kymo_Cnd=double(imread(condensin_name));end;     
    kymo_DNA=double(imread(dna_name));    
else                                %make two kymographs
    dna_pth=char(strcat(generaldatapth, Exp,'\', Channel_list(1)));
    condensin_pth=char(strcat(generaldatapth, Exp,'\', Channel_list(2)));    
    kymo_DNA=Build_kymo(dna_pth);
    kymo_Cnd=Build_kymo(condensin_pth);
end

% 2 do peak analysis on each profile of the kymographs


%loop analysis
    info_loop=Analyze_Loop(kymo_DNA);
    
if do_condensin
    info_condensin=PeakFit_kymo(kymo_Cnd,'flatbottom',2.7,initval);
    info_condensin=Label_Condensin_Loop(info_condensin,info_loop,initval);
    
end

%% save data   
    save(strcat(SaveName,   '_allresults.mat'),... 
                            'kymo_DNA',...
                            'info_loop');
    if do_condensin
        save(strcat(SaveName, '_allresults.mat'),... 
                              'kymo_Cnd',...
                              'info_condensin', '-append');
    end   
end
end


%% plot loop panel 1
if actions.plotmodus>0;    
for ee=1:1:LE
     
    close all;
    Exp=strcat('ROI',num2str(AllExp(ee)));    
    LoadName=char(strcat(outpath, Exp)); 
    load(strcat(LoadName, '_allresults.mat'));
    dna_name=char(strcat(generaldatapth, Exp, '\', Kymo_list,Dna_Kymo));
    condensin_name=char(strcat(generaldatapth, Exp,'\', Kymo_list,Condensin_Kymo));
    if exist(condensin_name)==2, do_condensin=1;, else  do_condensin=0;end  
    disp(strcat('Building&saving plots: Exps to work through:',num2str(LE-ee)));
    
     %% plot panel 1: overview
     if actions.plotmodus==1
     set(figure(1), 'visible', 'off');
    [rrd,ccd]=size(kymo_DNA);
     if do_condensin 
        subplot(2,3,1); 
            pcolor(kymo_Cnd); shading flat, colormap hot;
            title('Condensin'); ylabel('frame no.');              
        subplot(2,3,4);         
            pcolor(kymo_DNA); shading flat, colormap hot;
        title('DNA'); ylabel('frame no.');
     else
         subplot(1,3,1); 
            pcolor(kymo_DNA); shading flat, colormap hot;
        title('DNA'); ylabel('frame no.');
     end

    subplot(1,3,2);
        pairx=[info_loop.pos_loop_left' info_loop.pos_loop_right'];
        pairx2=[info_loop.pos_tether_left' info_loop.pos_tether_right'];
        pairy=[info_loop.frameno' info_loop.frameno'];      
        if do_condensin,
            sel=find(info_condensin.label_loopassociated==1);
            plot(info_condensin.pos_X_subpix(sel),info_condensin.pos_frameno(sel),'ro', 'MarkerSize',3); hold on; 
        end
        plot(info_loop.pos_loop_main,info_loop.frameno,'go'); hold on;
        plot(pairx',pairy','b-');
        plot(pairx2',pairy','k+');
        if do_condensin,legend('loop_condensins', 'loop main peak','loop edges','tether edges'); else
        legend('loop main peak','loop edges','tether edges'); end
        xlim([1 ccd]); ylim([1 rrd]);
        plot(pairx,pairy,'bo','MarkerSize',3,'MarkerFaceColor','b');    
        plot(info_loop.pos_loop_main,info_loop.frameno,'go','MarkerSize',3);    
        if do_condensin
            sel=find(info_condensin.label_loopassociated==1);
            plot(info_condensin.pos_X_subpix(sel),info_condensin.pos_frameno(sel),'ro', 'MarkerSize',3); hold on; 
        end   
        xlabel('position, pixels');
        ylabel('frame no.');
        hold off;
    
    subplot(1,3,3);
        plot(info_loop.frameno, info_loop.cont_loop_excess_raw, 'k-'); hold on;
        plot(info_loop.frameno, info_loop.cont_loop_excess, 'r-'); hold on;
        plot(info_loop.frameno, info_loop.cont_loop_mushroom, 'y-'); hold on;
        plot(info_loop.frameno, info_loop.cont_tether_left, 'm-'); hold on;
        plot(info_loop.frameno, info_loop.cont_tether_right, 'b-'); hold off;
        legend('loop only-raw','loop only', 'mid mushroom', 'left', 'right','Location','SouthOutside');
        xlabel('frame no.');
        ylabel('loopcontent, %');
        xlim([1 info_loop.frameno(end)]); 
        hold off;
  saveas(gcf,strcat(LoadName, '_plots.jpg'),'jpg'); 
  end
    
    %% plot panel 2: condensin label result 
 if do_condensin&&(actions.plotmodus==2)
        set(figure(2), 'visible', 'off');
  
    %         info_condensin
                    % pos_frameno
                    % pos_X_pix
                    % pos_X_subpix
                    % content_peakvals
                    % content_perspot_est
                    % content_perspot_meas
                    % label_OKspot
                    % label_looplabel
                    % label_generaledgelabel
    LC=length(info_condensin.pos_frameno);
    sel1=find((info_condensin.label_OKspot)==1); LOK=length(sel1);
    sel2=find(((info_condensin.label_OKspot)==1)...
              &(info_condensin.label_loopassociated==1)); LLP=length(sel2);
    sel3=find(((info_condensin.label_OKspot)==1)...
              &(info_condensin.label_generaledgelabel==3)); LLM=length(sel3);
    sel4=find(((info_condensin.label_OKspot)==1)...
              &((info_condensin.label_generaledgelabel==2)...
                |(info_condensin.label_generaledgelabel==4))); LLE=length(sel4);
     sel5=find(((info_condensin.label_OKspot)==1)...
              &((info_condensin.label_generaledgelabel==1)...
                |(info_condensin.label_generaledgelabel==5))); TTE=length(sel5);        
    sel6=find(((info_condensin.label_OKspot)==1)...
              &(info_condensin.label_generaledgelabel==0)); NN=length(sel6);                
    
          
          
     SpotSummaryLegend=([{'1.strengthOK'},{'2.loop'},...
                        {'3.loop edges'},{'4.loop brightest'} ,{'5.tether edges'},...
                        {'6.none'}]);   
                    clr=[{'k'},{'b'},{'r'},{'m'},{'y'},{'g'}];
    SpotSummaryVals=[LOK LLP  LLE LLM TTE NN];
    LS=length(SpotSummaryVals);
    subplot(1,3,1);
    axz=1:LS;
    for ii=1:LS
        bar(axz(ii),SpotSummaryVals(ii),char(clr(ii))); hold on;       
    end
    axis tight;
    legend(SpotSummaryLegend,'Location','NorthOutside');
    title(strcat('Spot categoriesof: ',Exp));
    binax=0:25:1000;
    
    subplot(3,2,2);
    Freelabels=info_condensin.content_perspot_est(sel6);
    hist_free=hist(Freelabels,binax);
    bar(binax,hist_free,'k'); xlim([0 1000]); 
    legend('free');
    
    subplot(3,2,4);
    Edgelabels=info_condensin.content_perspot_est(sel4);
    hist_edge=hist(Edgelabels,binax);
    bar(binax,hist_edge,'r'); xlim([0 1000]); 
    legend('loop edges');
    
     subplot(3,2,6);
    Mainlabels=info_condensin.content_perspot_est(sel3);
    hist_main=hist(Mainlabels,binax);
    bar(binax,hist_main,'m'); xlim([0 1000]); 
    legend('loop brightest');
    
    saveas(gcf,strcat(LoadName, '_spothistograms.jpg'),'jpg'); 
 end
end
end
end

function kymo=Build_kymo(pth)
    codepth=pwd;
    cd(pth); imlist=dir('*.tif'); cd(codepth);
    FramesNo=length(imlist);
    for jj=1:FramesNo  %for all frames        
        if mod(jj,50)==0, disp(strcat('Stack frames to load:',num2str(FramesNo-jj)));end
        source=strcat(pth,imlist(jj).name);
        pic=(double(imread(source)));
        if jj==1
             [rr,cc]=size(double(pic));
             kymo=zeros(FramesNo,rr);         
        end
        prf=sum(pic');
        prf=prf-min(prf);
        kymo(jj,:)=prf;   
    end
    
    
function condensin_info=PeakFit_kymo(kymo,peakfitoption, sigma,initval);
posses=[];
condensin_info=struct('pos_frameno',[]);
[FramesNo,~]=size(kymo);
cleankymo=kymo;
    switch peakfitoption
        case 'peeling2'
            for jj=1:FramesNo 
                prf=kymo(jj,:);
                prf=smooth(prf,4);
                %[peakprops,buildprf]=peel_peaks_from_profile(prf',2.7,1);
                [peakprops,buildprf,~, ~]=peel_peaks_from_profile_plusclusters(prf',sigma,1,'NonPeriodic');
                cleankymo(jj,:)=buildprf-median(buildprf);
                [betterpeaks, betterpeaksvals]= Refine_Peaks(prf,peakprops(:,3), 0);
                [LL,~]=size(peakprops);    
                posses=[posses; [jj+zeros(LL,1) betterpeaks peakprops(:,3)]];
            end
        case 'flatbottom'
             for jj=1:FramesNo 
                prf=kymo(jj,:); LP=length(prf); 
                 prf=JKD1_PRF_smooth(prf',4);
                 prf=prf-min(prf);
                [firstpeaks,betterpeaks, betterpeaksvals]=...
                 JKD1_PRF_get1Dpeaksflatbottom(prf,sigma,1,0);%data,sigs,refine,plotit
                dum=1;
                content=betterpeaksvals*((2*pi)^0.5*initval.psf_est);
                LL=length(betterpeaks); cont_m=zeros(LL,1); 
                if LP>0
                idxes=1:LP; 
                for pp=1:LL
                    idx=firstpeaks(pp);
                    sel=find((idxes>idx-2*initval.psf_est)&(idxes<idx+2*initval.psf_est));
                    cont_m(pp)=sum(prf(sel));
                    dum=1;
                end
                
                end
                
                posses=[posses; [jj+zeros(LL,1) betterpeaks firstpeaks, betterpeaksvals content cont_m]];
             end
        case 'edgedetect'        
            [stepfitim_c,stepfitim_r,stepcoords_r,stepcoords_c]=StepfindImage(kymo);
            stepcoords_c=Merge_Edges(stepcoords_c,stepfitim_c, kymo);  
            posses=[stepcoords_c(:,2) realshrinkfactor*stepcoords_c(:,3) realshrinkfactor*stepcoords_c(:,3)];
    end
    
    condensin_info.pos_frameno=posses(:,1)';
    condensin_info.pos_X_pix=posses(:,3)';
    condensin_info.pos_X_subpix=posses(:,2)';
    condensin_info.content_peakvals=posses(:,4)';
    condensin_info.content_perspot_est=posses(:,5)';
    condensin_info.content_perspot_meas=posses(:,6)';
    dum=1;
    
    
     function [kymo_shr, realshrinkfactor]=shrink_kymo(kymo, shrinkfactor);
     %shorten along x with about the psf in pixels
     
     [rr,cc0]=size(kymo);
     cc1=round(cc0/shrinkfactor);
     kymo_shr=zeros(rr,cc1);
     ax0=1:cc0;
     ax1=linspace(1,cc0,cc1);
     for ii=1:rr
         prf0=kymo(ii,:);
         prf1=interp1(ax0,prf0,ax1,'spline');
         kymo_shr(ii,:)=prf1;
     end
     realshrinkfactor=cc0/cc1;   
         
    
   function loopinfo=Analyze_Loop(kymo_DNA);
       %assume one dna loop ('cluster')
       % 1. set a left and right cluster border, based on treshold from
       % buildcurve
       %2. look for all spots (irrespective of cluster) that fall in this
       %slot. find leftmost and rightmost spot. The are the loop start and
       %stop positons, psf-corrected
       %get out:  
        
        %step 1: identify loop 
        kymo_DNA=pre_cook_kymo_DNA(kymo_DNA);
        fluo=Get_FluoLevelProps(kymo_DNA);
        
        %define a kymograph containing as much as possible only loop signal
        loop_kymo=kymo_DNA-fluo.level_looptreshold;
        loop_kymo(loop_kymo<0)=0;   %shave off level
        looplevel=mean(loop_kymo(loop_kymo~=0));       
        residu_kymo=kymo_DNA-loop_kymo;      
         
        %define a kymograph with fluorescent content only (no dark noise)
        content_kymo=kymo_DNA-fluo.level_darktreshold;
        content_kymo(content_kymo<0)=0;   %shave off level     
        [FramesNo,~]=size(loop_kymo);
        skips=1;
        for jj=1:skips:FramesNo  
            %loop analysis: get profiles
                prf_ori=kymo_DNA(jj,:);
                prf_ori=prf_ori-min(prf_ori);
            
                prf_loop=loop_kymo(jj,:); 
                prf_loop=prf_loop-min(prf_loop);
                
                prf_res=residu_kymo(jj,:);
                prf_res=prf_res-min(prf_res);
                
                prf_cont=content_kymo(jj,:);
                prf_cont=prf_cont-min(prf_cont);
      
    
            %get edges of tether
           [tetherstart,tetherstop,tetherok]=Get_edgeslength(prf_res,'tether');
           loopinfo.pos_tether_left(jj)=tetherstart;
           loopinfo.pos_tether_right(jj)=tetherstop;
            
     %% get loop structure
           %see if there is a loop too begin with and not wrongly (too
           %wide) detected
           loopinfo.frameno(jj)=jj;
           loopinfo.cont_loop_excess_raw(jj)=100*sum(prf_loop)/sum(prf_cont);
           if loopinfo.cont_loop_excess_raw(jj)>13    
                prf_loop=(smooth(prf_loop,4));            
                prf_loop=prf_loop';
                [lft_edge,rgt_edge,loopok]=Get_edgeslength(prf_loop,'loop');           
                loopinfo.loopok(jj)=loopok;
                if (lft_edge-rgt_edge)/(tetherstart-tetherstop) >0.7
                    loopinfo.loopok(jj)=0;
                end
           else
               loopinfo.loopok(jj)=0;
           end
           dum=1; 
 
            if loopinfo.loopok(jj)
                [peakprops,~]=peel_peaks_from_profile(prf_loop,2.7,0); %1 split profile into spot-components
                %[count PeakVal Xpos  Psf ThisSpotFraction CoveredFraction RelChange]];
                
                %2 find loop-encompassed component positions
                inpeaks_ix=find((peakprops(:,3)>lft_edge)&(peakprops(:,3)<rgt_edge));
                
                
                if ~isempty(inpeaks_ix);
                    inpeaks_posses=sort(peakprops(inpeaks_ix,3));
                    inpeakpos_leftmost=inpeaks_posses(1);
                    inpeakpos_righmost=inpeaks_posses(end);
                    inpeaks=peakprops(inpeaks_ix,:);
                    [~,idx_mx]=nanmax(inpeaks(:,5)); 
                    inpeaks_content=nansum(peakprops(inpeaks_ix,5));
                    
                    loopinfo.pos_loop_main(jj)=inpeaks(idx_mx,3);  %pos main blob                                   
                    loopinfo.pos_loop_left(jj)=inpeakpos_leftmost;          %pos left        
                    loopinfo.pos_loop_right(jj)=inpeakpos_righmost;         %pos left       
                    loopinfo.pos_loop_length(jj)=inpeakpos_righmost-inpeakpos_leftmost;
                    
                    %collect contents, all in percentages

                    loopinfo.cont_loop_excess(jj)=100*inpeaks_content*sum(prf_loop)/(sum(prf_loop+prf_res));                    
                    loopinfo.cont_tether_mid(jj)=100*sum(prf_cont(lft_edge+1:rgt_edge-1))/sum(prf_cont);;   %same, in and under loop
                    loopinfo.cont_loop_mushroom(jj)=loopinfo.cont_loop_excess(jj)+loopinfo.cont_tether_mid(jj);
                    loopinfo.cont_tether_left(jj)=100*sum(prf_cont(1:lft_edge))/sum(prf_cont); %all fluorescence left of loop                    
                    loopinfo.cont_tether_right(jj)=100*sum(prf_cont(rgt_edge:end))/sum(prf_cont); ; %same, right side                                    
                    loopinfo.cont_tether_leftonly(jj)=loopinfo.cont_tether_left(jj)+loopinfo.cont_tether_mid(jj)/2;
                    loopinfo.cont_tether_rightonly(jj)=loopinfo.cont_tether_right(jj)+loopinfo.cont_tether_mid(jj)/2;
                else
                    loopinfo=NaN_Loopinfo(loopinfo,jj);                   
                end
            else
                loopinfo=NaN_Loopinfo(loopinfo,jj);
            end
            loopinfo=orderfields(loopinfo);
            
            if 0
                %subplot(1,2,1);
                plot(prf_ori,'k-'); hold on;
                plot(prf_loop); hold on;                               
                plot(prf_cont,'r-');
                plot(prf_res,'b-');
                xlabel('positition, pixel units');
                ylabel('fluorescence intensity, a.u.');
                legend('\fontsize{8} ori','\fontsize{8} loop',...
                       '\fontsize{8} content','\fontsize{8} residu'); 
                %[~]=ginput(1);
                hold off;
                 subplot(1,2,2);
                 plot(loopinfo.frameno(1:skips:jj), loopinfo.cont_loop_excess_raw(1:skips:jj),'k');
                pause(0.01);
            end
            
            dum=1;
        end
        %[~]=ginput(1);
        
        
 function fluo=Get_FluoLevelProps(kymo);         
     kymo=double(kymo');  %position is vertical
    [rr,cc]=size(kymo);
    %1 First, we make estimates on intensity and noise of the background (outside
    %the DNA).'local background' is defined as the average of the outer two
    %image lines'. Then estimate dark (camera) noise via the spread in the difference between
    %neighbouring pixels
        fluo.level_dark=(mean(mean(kymo(1:2,:)))+mean(mean(kymo(rr-1:rr,:))))/2; 
        diftop=kymo(1:2,2:end)-kymo(1:2,1:end-1); 
        difbot=kymo(rr-2:rr,2:end)-kymo(rr-2:rr,1:end-1);
        dif=[diftop(:);  difbot(:)];
        fluo.noise_dark=std(dif)/2^0.5;
    
       %define an' surely inside tether' area
       
       tetherarea=kymo(20:end-20,:);
       fluo.tetherlevel=median(tetherarea(:));
    
    %define as 'fluorescence' those pixels sufficiently above the darklevel.
    %Note this may not be representative for the outline of the bacterium,
    %since there is some blurring and we want to measure all fluorescence
    fluo.level_darktreshold=fluo.level_dark+2*fluo.noise_dark;
    fluo.level_tethertreshold=fluo.tetherlevel-2*fluo.noise_dark;
    fluo.level_looptreshold=fluo.tetherlevel+2*fluo.noise_dark;
    
   dum=1;

   
 function [curvestart,curvestop,curveok]=Get_edgeslength(prf_res, type_of_profile);
    %function uses 'shaved off' profile to find start and stop
    tresval=0.4;
    Lp=length(prf_res); axz=1:Lp;
        
    [lo_L,ixL]=min(prf_res(1:ceil(Lp/2)));
    [lo_R,ixR]=min(prf_res(ceil(Lp/2):end)); ixR=ixR+ceil(Lp/2)-1;
    slopefit=polyval(polyfit([ixL ixR],[lo_L lo_R],1),axz);
    
    prf_res_ft=prf_res-slopefit;
    
    switch type_of_profile
        case 'tether' %main level, assuming most of profile is tether 
            main_lev=median(prf_res_ft);
            %subplot(2,2,2);
        case 'loop' %%FWHM level, assuming relatively compact loop 
            main_lev=0.5*max(prf_res_ft);
            %subplot(2,2,4);
    end
    lo=min(prf_res_ft);
    sel=find(prf_res_ft>(lo+tresval*(main_lev-lo)));
    
    if ~isempty(sel);
        curvestart=min(sel);     
        curvestop=max(sel);
        curveok=1;
    else
        curvestart=NaN;     
        curvestop=NaN;
        curveok=0;
    end

    if 0       
        plot(prf_res,'b-'); hold on;
        plot(slopefit,'b-');
        plot(prf_res_ft,'r');       
        plot(0*prf_res+main_lev,'r-');
        stem(curvestart,prf_res_ft(curvestart),'ro');
        stem(curvestop,prf_res_ft(curvestop),'ro');
        legend( '\fontsize{6} residu','\fontsize{6} minima fit',...
                '\fontsize{6} corrected','\fontsize{6} mean',...
                '\fontsize{6}edges');
        xlabel('positition, pixel units');
        ylabel('fluorescence intensity, a.u.');
        pause(0.5);        
        [~]=ginput(1);
        hold off;
    end
    
    
        
  function kymo_out=pre_cook_kymo_DNA(kymo);
      %clean up: background slope, bleaching
      %1 flatten background left-right
      
      prf_med=median(kymo);
      allminval=min(prf_med);
      [ff,Lp]=size(kymo);
      xax=1:Lp;
        
      [lo_L,ixL]=min(prf_med(1:ceil(Lp/2)));
      [lo_R,ixR]=min(prf_med(ceil(Lp/2):end)); ixR=ixR+ceil(Lp/2)-1;
      slopeline=polyval(polyfit([ixL ixR],[lo_L lo_R],1),xax);
      slopeplane=repmat(slopeline,ff,1);
      kymo=kymo-slopeplane+allminval;
      
      
      %2 crude bleach correct
      tax=1:ff;
      bleachline=mean([kymo(1:ff,1) kymo(1:ff,end)]');
      bleachline_nrm=bleachline/mean(bleachline);
      bleachslope=polyval(polyfit(tax,bleachline_nrm,1),tax);
      bleachplane=repmat(bleachslope',1,Lp);
      
      kymo=kymo./bleachplane;
      kymo_out=JKD2_IM_smoothJK(kymo,3);
      
      for ii=1:ff
          kymo_out(ii,:)=kymo_out(ii,:)-min(kymo_out(ii,:));
      end
      
      
      if 0          
        subplot(1,2,1);
            plot(tax(2:end),bleachline_nrm(2:end)); hold on;
            plot(tax(2:end),bleachslope(2:end),'r');
            title('bleach correction');
            legend('mean edges, norm.','linear fit');
            xlabel('time, frames');
            ylabel('norm.intensity, a.u.');
            xlim([1 ff]);
        subplot(1,2,2);           
            plot(prf_med,'b-'); hold on;
            plot(slopeline,'b-');
            plot(prf_med-slopeline+allminval,'r'); 
            title('flattening');
            legend('median','minima fit','corrected');
            xlabel('position, pixel units');
            ylabel('fluorescence intensity, a.u.');
            pause(0.5);        
            [~]=ginput(1);
            hold off
      end
      
      function loopinfo=NaN_Loopinfo(loopinfo,jj);
        loopinfo.pos_loop_main(jj)=NaN;                                   
        loopinfo.pos_loop_left(jj)=NaN;       
        loopinfo.pos_loop_right(jj)=NaN;     
        loopinfo.pos_loop_length(jj)=NaN;
        loopinfo.cont_loop_excess(jj)=NaN;                    
        loopinfo.cont_tether_mid(jj)=NaN;
        loopinfo.cont_loop_midplusexcess(jj)=NaN;
        loopinfo.cont_tether_left(jj)=NaN;                 
        loopinfo.cont_tether_right(jj)=NaN;                                
        loopinfo.cont_tether_leftonly(jj)=NaN;
        loopinfo.cont_tether_rightonly(jj)=NaN;
        
        function info_condensin=Label_Condensin_Loop(info_condensin,info_loop,initval);
        %this fuction screens and labels the detected condensin spots
%         info_condensin
                    % pos_frameno
                    % pos_X_pix
                    % pos_X_subpix
                    % content_peakvals
                    % content_perspot_est
                    % content_perspot_meas
         info_condensin.label_OKspot=0*info_condensin.pos_frameno; 
         info_condensin.label_loopassociated=0*info_condensin.pos_frameno;
         info_condensin.label_generaledgelabel=0*info_condensin.pos_frameno;
         
         Ic=info_condensin.content_perspot_meas;
         %first, get some measures
        [flag,cleandata]=JKD1_PRF_outlier_flag(Ic,4,0.7,'all',0); 
        sel=find(flag==1);
        info_condensin.label_OKspot(sel)=1;
        
        %now, per time point
        FF=max(info_loop.frameno);
        for ii=1:FF
            fi_L=find(info_loop.frameno==ii);
            fi_C=find(info_condensin.pos_frameno==ii);
            
            if ~isempty(fi_L)&~isempty(fi_C);
            
                %loop positions of interest
                L_lft=info_loop.pos_loop_left(fi_L);
                L_mn=info_loop.pos_loop_main(fi_L);
                L_rgt=info_loop.pos_loop_right(fi_L);
                Th_lft=info_loop.pos_tether_left(fi_L);
                Th_rgt=info_loop.pos_tether_right(fi_L);
                Dna_geometry=[Th_lft L_lft L_mn L_rgt Th_rgt];
                
                Cxx=info_condensin.pos_X_subpix(fi_C); Lc=length(Cxx);
                for jj=1:Lc
                    oridx=fi_C(jj); %original indices of condensins;
                    Cx=Cxx(jj);
                    dd=abs(Dna_geometry-Cx);
                    [val,labl]=min(dd);  %nearest of these
                    if (Cx>L_lft-initval.psf_est)&(Cx<L_rgt+initval.psf_est)
                        info_condensin.label_loopassociated(oridx)=1;
                    end
                    if val<initval.psf_est;  %localized with
                        info_condensin.label_generaledgelabel(oridx)=labl;
                    end
                 dum=1;   
                end
                dum=1;
            end
        end
        dum=1;
        
        
        
            
            
            
        