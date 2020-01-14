function A020_Condensin_and_plectonemes_get_kymographs_and_positions(init,expi);
%JWJK_A:-------------------------------------------------------------------
%Summary: %This function builds kymographs from images 
%and analyzes them for peaks associated with DNA plectonemes and condensin
%JacobKers2019
%Approach: kymographs are made from movie stacks; peak detection is 
%performed; the positions of these peaks (condnesin and plectonemes) are
%related to each other
%run actions can be set or unset, to save time.
%Input: exp info and movie stacks.
%Output: kymographs and spot information
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_A-------------------------------------------------------------------
close all;


%% 3 choose the run actions; 
%each step stores new data&reloads from former step 
actions.buildkymographs=1;  %make raw kymographs
actions.peakdetection=1;    %detect peaks; convert to genomic percentage and condensin counts
actions.plot=1;   %plectoneme/condensin relations etc
    

%% 1) Set common paths; use standardized naming
datapathin=init.datapathin;
datapathout=init.datapathout;
expname=init.expname;
AllExp=init.AllExp;

generaldatapth=[datapathin,expname,'\'];
outpath=strcat(datapathout, 'matlabresults\',expname,'\');
if ~isdir(outpath), mkdir(outpath); end



psf_est=init.psf_est;
LE=length(AllExp);  %for all experiments
for roi=1:LE  
Exp=strcat(init.roidirname,num2str(AllExp(roi)));
    expinfo=A002_Condensin_with_plectonemes_expinfo(expi,AllExp(roi));
    SaveName=char(strcat(outpath,'EKMcp_A020_',Exp)); 
    
%% Kymographs
if actions.buildkymographs     
    if mod(roi,1)==0, disp(strcat('Building kymograph:',expname,':Exps to work through:',num2str(LE-roi+1)));end     
    %% 1)make two kymographs 
        Channel_list=[{'DNA\'}, {expinfo.labelname}];  
        dna_pth=char(strcat(generaldatapth, Exp,'\', Channel_list(1)));
        condensin_pth=char(strcat(generaldatapth, Exp,'\', Channel_list(2)));
        
        kymo_DNA=kym_build_kymo_from_movie(dna_pth,expinfo);       
        kymo_Cnd=kym_build_kymo_from_movie(condensin_pth,expinfo);             
        save(strcat(SaveName, '_allresults.mat'), 'kymo_DNA','kymo_Cnd');
        
        if 1 %optional back-saving a kymograph  for later use (for example, using 'boxtrack')
          dlmwrite([generaldatapth, Exp,'\kymograph\EKMcp_A020_Kymograph_MukBEF.txt'],kymo_Cnd);
        end
        
        if 1
        figure(64); 
            subplot(1,2,1); pcolor(kymo_DNA); shading flat; colormap hot; 
            title('DNA');
            subplot(1,2,2); pcolor(kymo_Cnd); shading flat; colormap hot; 
            title('Condensin');
            %[~]=ginput(1); 
            saveas(gcf,strcat(SaveName, '_kymographs.jpg'),'jpg'); 
            pause(2);
            close(gcf);
        end
    end   

%% Analysis   
if actions.peakdetection        
    if mod(roi,1)==0, disp(strcat('Analyzing:',expname,':Exps to work through:',num2str(LE-roi+1)));end 
    load(strcat(SaveName, '_allresults.mat'));
    %% 2 peak analysis on each profile of the kymographs
        kymo_DNA=kym_smooth_time_or_place(kymo_DNA,init.t_smooth,'time');
        kymo_Cnd=kym_smooth_time_or_place(kymo_Cnd,init.t_smooth,'time');  
        kymo_DNA=kym_smooth_time_or_place(kymo_DNA,init.x_smooth,'place');
        kymo_Cnd=kym_smooth_time_or_place(kymo_Cnd,init.x_smooth,'place');  
    
    
    %get some general properties  
        %Condensin
        levels_Cnd=kym_get_signal_levels(kymo_Cnd,init.tresholdsigmas);
        kymo_Cnd_peaks=kym_keep_only_peaks(kymo_Cnd,levels_Cnd); %shave off loops  
        info_Cnd=kym_peakfitperkymographline(kymo_Cnd,'flatbottom', psf_est,expinfo.tres_pk_Cnd); %build spot info!
        
        
        
        %DNA
        levels_DNA=kym_get_signal_levels(kymo_DNA,init.tresholdsigmas);  
        kymo_DNA=kym_convert_dna_kymo(kymo_DNA,levels_DNA);  %convert_DNA counts to_genomic_percentage;
        levels_DNA=kym_get_signal_levels(kymo_DNA,init.tresholdsigmas);  % repeat levels 
        kymo_DNA_peaks=kym_keep_only_peaks(kymo_DNA,levels_DNA); %shave off loops  
        info_DNA=kym_peakfitperkymographline(kymo_DNA_peaks,'just_treshold', psf_est,expinfo.tres_pk_DNA);   %build spot info!
        
        [info_DNA,info_Cnd]=kym_get_length_and_densities(info_DNA,info_Cnd,kymo_DNA);
     
 
        save(strcat(SaveName,'_allresults.mat'),'info_DNA','info_Cnd',...
                            'kymo_DNA_peaks', 'kymo_Cnd_peaks', '-append');    
end


%% plotting
if actions.plot>0       
    close all;
    Exp=strcat('ROI',num2str(AllExp(roi)));  
    disp(strcat('Building&saving plots: Exps to work through:',num2str(LE-roi)));
    load(strcat(SaveName, '_allresults.mat')); 
    
    
     %% plot panel 1: overview
     %set(figure(1), 'visible', 'off');
        plotpic_DNA=kymo_DNA_peaks;
        plotpic_DNA(plotpic_DNA==0)=NaN;
        plotpic_DNA=log10(plotpic_DNA);
        plotpic_Cnd=kymo_Cnd_peaks;
        plotpic_Cnd(plotpic_Cnd==0)=NaN;
        plotpic_Cnd=log10(plotpic_Cnd);
        
%         subplot(1,3,1); 
%        
%             pcolor(plotpic_DNA); shading flat; colormap jet; hold on;
%             title('DNA');
%         subplot(1,3,2); 
%             pcolor(plotpic_Cnd); shading flat; colormap jet; hold on;
%             title('Condensin');

       subplot(1,3,2); 
            plot(info_DNA.pos_X_subpix+expinfo.channelshift, info_DNA.pos_frameno, 'bo','Markersize',4); hold on;
            plot(info_Cnd.pos_X_subpix, info_Cnd.pos_frameno, 'ro','Markersize',2);
            xlim([0 length(kymo_DNA(1,:))]);
            ylim([0 length(kymo_DNA(:,1))]);
            title(Exp);
            legend('DNA','Condensin')
            %[sh,~,~]=ginput(2); 
            %xshift=sh(2)-sh(1) %output shift
            saveas(gcf,strcat(SaveName, '_plectonemecounts.jpg'),'jpg');    
            pause(1);
            close(gcf); 
            
            if 1
        figure(165); 
            subplot(1,2,1); pcolor(kymo_DNA); shading flat; colormap hot; hold on;
            plot(info_DNA.pos_X_subpix+expinfo.channelshift, info_DNA.pos_frameno, 'wo','Markersize',1); hold on;           
            xlim([0 length(kymo_DNA(1,:))]);
            ylim([0 length(kymo_DNA(:,1))]);
            title('DNA');
            subplot(1,2,2); pcolor(kymo_Cnd); shading flat; colormap hot; hold on;
            plot(info_Cnd.pos_X_subpix, info_Cnd.pos_frameno, 'wo','Markersize',1);
            xlim([0 length(kymo_DNA(1,:))]);
            ylim([0 length(kymo_DNA(:,1))]);
            title('Condensin');
            pause(0.1);
            %[~]=ginput(1); 
            saveas(gcf,strcat(SaveName, '_kymographs_overlays.jpg'),'jpg'); 
            
            close(gcf);
        end
end
end


  



        