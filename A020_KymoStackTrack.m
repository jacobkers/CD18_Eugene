
function A020_KymoStackTrack
%JWJK_A:-------------------------------------------------------------------
%Title:A020_KymoStackTrack
%
%Summary: %This function builds a kymograph from images 
%and analyzes it
%JacobKers2018
%
%:JWJK_A-------------------------------------------------------------------
actions.analyze=1;


close all;
loadImageJ_kymograph=1;
datapath='D:\jkerssemakers\_Data\CD\2018_Eugene\';
exprun=1;
switch exprun
    case 1
        generaldatapth=[datapath,'2018_08_01 Pilot Runs\'];
        AllExp=[3 1];  %numbers of various rois   
    case 2
        generaldatapth=[datapath,'\2018_09_24 More_molecules\'];
        AllExp=[2 4 5 6 7 8 9 10 11 12 13 14 15 16];  %paths to various rois  
        %AllExp=[7 8 9 10 11 12 13 14 15 16];  %paths to various rois   
        AllExp=[13];  %paths to various rois   
end

outpath=strcat(generaldatapth, 'matlabresults\');
 

if ~isdir(outpath), mkdir(outpath); end

%% standardized subdirectory names 
Channel_list=[{'DNA\'}, {'Condensin\'}];     %The two subchannels
Kymo_list='kymo_ImageJ\'; %ImageJ-made kymographs
Condensin_Kymo='Kymograph_Condensin.tif';           %if you use it
Dna_Kymo='Kymograph_DNA.tif';                %if you use it';


%% main loop
LE=length(AllExp);  %for all experiments
if actions.analyze
for ee=1:1:LE
if mod(ee,1)==0, disp(strcat('Analyzing: Exps to work through:',num2str(LE-ee)));end 
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
if do_condensin
    posses_Cnd=PeakFit_kymo(kymo_Cnd,'flatbottom',4);
end

%2a do loop analysis
loopinfo=Analyze_Loop(kymo_DNA);

%% save data   
    save(strcat(SaveName,   '_allresults.mat'),... 
                            'kymo_DNA',...
                            'loopinfo');
    if do_condensin
        save(strcat(SaveName, '_allresults.mat'),... 
                              'kymo_Cnd',...
                              'posses_Cnd', '-append');
    end   
end
end

%% plot loop

for ee=1:1:LE
    Exp=strcat('ROI',num2str(AllExp(ee)));    
    LoadName=char(strcat(outpath, Exp)); 
    disp(strcat('Plotting: Exps to work through:',num2str(LE-ee)));
    load(strcat(LoadName, '_allresults.mat'));
    dna_name=char(strcat(generaldatapth, Exp, '\', Kymo_list,Dna_Kymo));
    condensin_name=char(strcat(generaldatapth, Exp,'\', Kymo_list,Condensin_Kymo));
    if exist(condensin_name)==2, do_condensin=1;, else  do_condensin=0;end   
    
    
    figure(1);
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
        pairx=[loopinfo.left' loopinfo.right'];
        pairx2=[loopinfo.tetherleft' loopinfo.tetherright'];
        pairy=[loopinfo.fr' loopinfo.fr'];      
        if do_condensin,plot(posses_Cnd(:,2),posses_Cnd(:,1),'ro', 'MarkerSize',6); hold on; end
        plot(loopinfo.main,loopinfo.fr,'go'); hold on;
        plot(pairx',pairy','b-');
        plot(pairx2',pairy','k+');
        if do_condensin,legend('condensin', 'loop main peak','loop edges','tether edges'); else
        legend('loop main peak','loop edges','tether edges'); end
        xlim([1 ccd]); ylim([1 rrd]);
        plot(pairx,pairy,'bo','MarkerSize',3,'MarkerFaceColor','b');    
        plot(loopinfo.main,loopinfo.fr,'go','MarkerSize',3);    
        if do_condensin,
            plot(posses_Cnd(:,2),posses_Cnd(:,1),'ro', 'MarkerSize',6); hold on; 
        end   
        xlabel('position, pixels');
        ylabel('frame no.');
        hold off;
    
    subplot(1,3,3);
        plot(loopinfo.fr, loopinfo.perc, 'r-'); hold on;
        plot(loopinfo.fr, loopinfo.perc_mid, 'm-'); hold on;
        plot(loopinfo.fr, loopinfo.perc_left, 'k-'); hold on;
        plot(loopinfo.fr, loopinfo.perc_right, 'b-'); hold off;
        legend('loop only', 'loop section', 'left section', 'right section','Location','SouthOutside');
        xlabel('frame no.');
        ylabel('loopcontent, %');
        xlim([1 loopinfo.fr(end)]); 
        hold off;
    saveas(gcf,strcat(LoadName, '_plots.jpg'),'jpg');      
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
    
    
function posses=PeakFit_kymo(kymo,peakfitoption, sigma);
posses=[];
[FramesNo,~]=size(kymo);
cleankymo=kymo;
    switch peakfitoption
        case 'peeling'
            for jj=1:FramesNo 
                prf=kymo(jj,:);
                prf=smooth(prf,4);
                %[peakprops,buildprf]=peel_peaks_from_profile(prf',2.7,1);
                [peakprops,buildprf,clusterprops, allclusters]=peel_peaks_from_profile_plusclusters(prf',2.7,1,'NonPeriodic');
                cleankymo(jj,:)=buildprf-median(buildprf);
                [betterpeaks, betterpeaksvals]= Refine_Peaks(prf,peakprops(:,3), 0);
                [LL,~]=size(peakprops);    
                posses=[posses; [jj+zeros(LL,1) betterpeaks peakprops(:,3)]];
            end
        case 'flatbottom'
             for jj=1:FramesNo 
                prf=kymo(jj,:);
                 prf=JKD1_PRF_smooth(prf',4);
                [firstpeaks,betterpeaks, betterpeaksvals]=...
                 JKD1_PRF_get1Dpeaksflatbottom(prf,sigma,1,0);%data,sigs,refine,plotit
                dum=1;
                LL=length(betterpeaks);
                posses=[posses; [jj+zeros(LL,1) betterpeaks firstpeaks]];
             end
        case 'edgedetect'        
            [stepfitim_c,stepfitim_r,stepcoords_r,stepcoords_c]=StepfindImage(kymo);
            stepcoords_c=Merge_Edges(stepcoords_c,stepfitim_c, kymo);
            
            posses=[stepcoords_c(:,2) realshrinkfactor*stepcoords_c(:,3) realshrinkfactor*stepcoords_c(:,3)];
    end
    
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
        for jj=1:FramesNo 
            
            
            %loop analysis: get profiles
            prf_ori=kymo_DNA(jj,:);
            if 1
                prf_loop=loop_kymo(jj,:);
                prf_res=residu_kymo(jj,:);
                prf_cont=content_kymo(jj,:);
            else
                %[prf_loop,prf_res,prf_cont]=Split_profile(prf_ori);
            end       
            if 0
                plot(prf_ori,'k-'); hold on;
                plot(prf_loop); hold on;                               
                plot(prf_cont,'r-');
                plot(prf_res,'b-');
                xlabel('positition, pixel units');
                ylabel('fluorescence intensity, a.u.');
                legend('ori','loop','content','residu');               
                [~]=ginput(1);
                hold off
            end
            
           [tetherstart,tetherstop]=Get_tetherlength(prf_res);
           loopinfo.tetherleft(jj)=tetherstart;
           loopinfo.tetherright(jj)=tetherstop;
            
            %get loop structure
            prf_loop=(smooth(prf_loop,4));
            [peakprops,buildcurve]=peel_peaks_from_profile(prf_loop',2.7,0);
            prf_loop=prf_loop';
            % 1. set a left and right cluster border
            tresh=0.25*looplevel;
            sel=find(prf_loop>tresh);
            
            loopinfo.fr(jj)=jj;              
            if ~isempty(sel)
                %get the left and right edge of the loop              
                lft_edge=sel(1);
                rgt_edge=sel(end);
                               
                %2 find encompassed spot positions
                inpeaks_ix=find((peakprops(:,3)>lft_edge)&(peakprops(:,3)<rgt_edge));
                if ~isempty(inpeaks_ix);
                    inpeaks=peakprops(inpeaks_ix,:);
                    [~,idx_mx]=nanmax(inpeaks(:,5));
                    mainpeakx=inpeaks(idx_mx,3);               
                    loopinfo.main(jj)=mainpeakx;                
                    [inpeaks_x,sortidx]=sort(inpeaks(:,3));                
                    loopinfo.left(jj)=inpeaks_x(1);
                    loopinfo.right(jj)=inpeaks_x(end);                
                    loopinfo.length(jj)=inpeaks_x(end)-inpeaks_x(1);
                    inpeaks_content=nansum(peakprops(inpeaks_ix,5));
                    loopinfo.perc(jj)=100*inpeaks_content*sum(prf_loop)/(sum(prf_loop+prf_res));                    
                    leftperc=100*sum(prf_cont(1:lft_edge))/sum(prf_cont);
                    midperc=100*sum(prf_cont(lft_edge:rgt_edge))/sum(prf_cont);
                    rightperc=100*sum(prf_cont(rgt_edge:end))/sum(prf_cont);                   
                    loopinfo.perc_left(jj)=leftperc;  %all fluorescence left of loop
                    loopinfo.perc_mid(jj)=midperc;   %same, in and under loop
                    loopinfo.perc_right(jj)=rightperc; %same, right side                   
                  else
                    loopinfo.main(jj)=NaN;
                    loopinfo.left(jj)=NaN;
                    loopinfo.right(jj)=NaN;
                    loopinfo.length(jj)=NaN;
                    loopinfo.perc(jj)=NaN;     
                    loopinfo.perc_left(jj)=NaN;  %all fluorescence left of loop
                    loopinfo.perc_mid(jj)=NaN;   %same, in and under loop
                    loopinfo.perc_right(jj)=NaN; %same, right side
                    
                end
            end
        end
        
        
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

   
 function [tetherstart,tetherstop]=Get_tetherlength(prf_res);
    %function uses 'shaved off' profile to find start and stop
    tresval=0.4;
    Lp=length(prf_res); axz=1:Lp;
        
    [lo_L,ixL]=min(prf_res(1:ceil(Lp/2)));
    [lo_R,ixR]=min(prf_res(ceil(Lp/2):end)); ixR=ixR+ceil(Lp/2)-1;
    slopefit=polyval(polyfit([ixL ixR],[lo_L lo_R],1),axz);
    
    prf_res_ft=prf_res-slopefit;
    mid_lev=mean(prf_res_ft); %assuming most of profile is tether
    
    lo=min(prf_res_ft);
    sel=find(prf_res_ft>(lo+tresval*(mid_lev-lo)));
    
    if ~isempty(sel);
        tetherstart=min(sel);     tetherstop=max(sel);
    else
        tetherstart=NaN;     tetherstop=NaN;
    end
        
    
    if 0       
        plot(prf_res,'b-'); hold on;
        plot(slopefit,'b-');
        plot(prf_res_ft,'r');       
        plot(0*prf_res+mid_lev,'r-');
        stem(tetherstart,prf_res_ft(tetherstart),'ro');
        stem(tetherstop,prf_res_ft(tetherstop),'ro');
        legend('residu','minima fit','corrected','mean','edges');
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