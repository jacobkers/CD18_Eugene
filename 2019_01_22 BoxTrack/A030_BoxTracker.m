function A030_BoxTracker(initval, expname)
%JWJK_A:----[add ABCorC*----------------------------------------------------
%Title: tracking a user-selected peak through a kymograph
%Summary: part of a shell program; it uses a list of user-clicked start-and
%stop coordinates to identify the motion of local maxima in a kymograph.
%Approach: tracking is performed between a start and stop time; peaks are 
%analyzed by position and content within a box of user-defined width. 
%Before the start time, a %motionless-content measurement is done to detect 
%low-peak intensity.  %Tracking alows for a 'look-ahead' averaging for increased 
%robustness of tracking.
%Input: list of start-and stop coordinated provided by the 'A28' clicker
%Output: a structure 'looptraces' containing content properties per peak
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-----[add ABCorC*---------------------------------------------------
close all;
addpath(genpath(pwd));

%% get the ROIs to analyze
Lroi=length(initval.roilist);
%%  work trhough all ROIs
for roii=1:Lroi
    close all;
    figure(roii);
    
    %% load a region-of-interest kymograph
    savit=1;
    roino=initval.roilist(roii);
    Exp=strcat('ROI',num2str(roino));
    LoadName=char(strcat(initval.expi_outpath, 'EKMbt_A028_',Exp,'_clickinfo'));
    load(LoadName);
    SaveName=char(strcat(initval.expi_outpath, 'EKMbt_A028_',Exp));
    
    %source kymograph
    datainpath=strcat(initval.expi_inpath,initval.roilabel, num2str(roistartstop.roino),initval.kymodir);       
    source=[datainpath, initval.kymofile];
    trackmap=dlmread(source);
   
    
    
    %% get general properties, such as tether edges
    [ff,cc]=size(trackmap);   
    [tetherstart,tetherstop]=kym_get_tetheredges(trackmap);
    initval.tetherlevel=kym_get_tetherlevel(trackmap,tetherstart,tetherstop);
    initval.tetherstart=tetherstart;
    initval.tetherstop=tetherstop;       
    Nloops=length(roistartstop.start_x);
    looptraces=struct('Lx',[]);
    
    %% initialize loop info
    if 1
    for jj=1:Nloops 
        xj=roistartstop.start_x(jj); %first x for trace;
        fj=roistartstop.pre_t(jj); %first frame;
        %get first box
        boxprops.boxhalfwidth=roistartstop.loopanalysishalfwidth(jj);
        boxprops.cutlevel=initval.tetherlevel;
        boxprops.lookahead=0;
        [~,Ij,~]=prf_get_box_intensity(trackmap,xj,fj,boxprops); %first I
        looptraces(jj).Lx=xj;
        looptraces(jj).frame=fj;
     
    end
    end
    %% now, analyze the kymograph
    for ii=1:ff
        if (mod(ii,300)==0), disp(num2str(ii));end;
        looptraces=kym_analyze_traces(looptraces,trackmap,roistartstop,ii,initval);
    end
    if savit
        save(strcat(SaveName, '_pairtrackresults.mat'),... 
                                  'looptraces','initval');
    end
    
    %Plot_menu_I(trackmap, looptraces,savit, SaveName,Nloops);
    Plot_menu_III(trackmap, looptraces,savit, SaveName,Nloops);
    Plot_menu_II(trackmap, looptraces,savit, SaveName,Nloops);
    
end



    


function Plot_menu_I(trackmap, looptraces,savit, SaveName,Nloops);
    close all;
    figure;
    subplot(2,1,1);
    pcolor(trackmap'); colormap bone, shading flat; hold on;
    colormap(bone); 
    for jj=1:Nloops 
        plot(looptraces(jj).frame+0.5,looptraces(jj).Lx+0.5,'r-');
        plot(looptraces(jj).frame+0.5,looptraces(jj).curvestart+0.5,'y-');
        plot(looptraces(jj).frame+0.5,looptraces(jj).curvestop+0.5,'y-');       
    end 
    title('traces')
    ylabel('frameno, a.u.');
    xlabel('position, a.u.');
     subplot(2,1,2);
    for jj=1:Nloops
        plot(looptraces(jj).frame,looptraces(jj).Zloop.I_mid,'y-'); hold on;
        plot(looptraces(jj).frame,looptraces(jj).Zloop.I_left,'m-'); hold on;
        plot(looptraces(jj).frame,looptraces(jj).Zloop.I_right,'b-'); hold on;
        plot(looptraces(jj).frame,looptraces(jj).Zloop.I_checksum,'k-'); hold on;
        xlim([1 length(trackmap(:,1))]);
    end
    title('Z-loop');
    legend('loop','left','right');
     xlabel('frameno, a.u.');
     ylabel('intensity, % of total');
      if savit
        saveas(gcf, strcat(SaveName, '_Z_loop_pairtrackresults.svg'));
        saveas(gcf, strcat(SaveName, '_Z_loop_pairtrackresults.jpg'));
        saveas(gcf, strcat(SaveName, '_Z_loop_pairtrackresults.fig'));
        pause(1);
        close(gcf);
      end
 
      
function Plot_menu_II(trackmap, looptraces,savit, SaveName,Nloops);
    close all;
    figure;
    subplot(1,2,1);
    pcolor(trackmap); colormap bone, shading flat; hold on;
    colormap(bone); 
    for jj=1:Nloops 
        plot(looptraces(jj).Lx+0.5,looptraces(jj).frame+0.5,'y');       
    end   
     title('traces')
        ylabel('frameno, a.u.');
        xlabel('position, a.u.');
    for jj=1:Nloops
        subplot(2,2,2*jj);
        plot(looptraces(jj).frame,looptraces(jj).Regloop.I_hat,'k-','LineWidth',1); hold on;
        plot(looptraces(jj).frame,looptraces(jj).Regloop.I_left_neighbour,'b-'); hold on;
        plot(looptraces(jj).frame,looptraces(jj).Regloop.I_right_neighbour,'r-'); hold on;
        title(['loop',num2str(jj),'Nb-left=bl,rgt=rd,lp=blck']);    
        xlim([1 length(trackmap(:,1))]);
        ylim([0 100]);
        xlabel('frameno, a.u.');
        ylabel('intensity, % of total');
    end
    %legend('hat','left-full','right-full');
     
     
    
           
     if savit
        saveas(gcf, strcat(SaveName, '_regular_loop_segments.svg'));
        saveas(gcf, strcat(SaveName, '_regular_loop_segments.jpg'));
        saveas(gcf, strcat(SaveName, '_regular_loop_segments.fig'));
        pause(1);
        close(gcf);
     end      
 
function Plot_menu_III(trackmap, looptraces,savit, SaveName,Nloops);
    close all;
    figure;
    subplot(1,2,1);
    pcolor(trackmap); colormap bone, shading flat; hold on;
    colormap(bone); 
    for jj=1:Nloops 
        plot(looptraces(jj).Lx+0.5,looptraces(jj).frame+0.5,'y');       
    end   
     title('traces')
        ylabel('frameno, a.u.');
        xlabel('position, a.u.');
        
    for jj=1:Nloops
        subplot(2,2,2*jj);
        plot(looptraces(jj).frame,looptraces(jj).Regloop.I_hat,'k-','LineWidth',1); hold on;
        plot(looptraces(jj).frame,looptraces(jj).Regloop.I_left_full,'b-'); hold on;
        plot(looptraces(jj).frame,looptraces(jj).Regloop.I_right_full,'r-'); hold on;
        title(['loop',num2str(jj),'full-left=bl,rgt=rd,lp=blck']); 
        xlim([1 length(trackmap(:,1))]);
        ylim([0 100]);
        xlabel('frameno, a.u.');
        ylabel('intensity, % of total');
    end
    %legend('hat','left-full','right-full');
     
     
    
           
     if savit
        saveas(gcf, strcat(SaveName, '_regular_loop_full.svg'));
        saveas(gcf, strcat(SaveName, '_regular_loop_full.jpg'));
        saveas(gcf, strcat(SaveName, '_regular_loop_full.fig'));
        pause(1);
        close(gcf);
     end