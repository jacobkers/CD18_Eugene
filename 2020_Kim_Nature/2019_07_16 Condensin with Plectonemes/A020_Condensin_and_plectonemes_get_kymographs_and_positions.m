function A020_Condensin_and_plectonemes_get_kymographs_and_positions(init,expi,usr,actions);
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
actions.backsaving=0;    %optional back-saving a kymograph  for later use; uses specific names! check the lines


%% 1) Set common paths; use standardized naming
datapathin=init.datapathin;
datapathout=init.datapathout;
expname=init.expname;
AllExp=init.AllExp;

generaldatapth=[datapathin,expname,'\'];
outpath0=strcat(datapathout, 'matlabresults\',init.expname,'\');
outpath1=strcat(datapathout, 'matlabresults\',init.expname,'\A020_kymographs\');
if ~isdir(outpath1), mkdir(outpath1); end

psf_est=init.psf_est;
LE=length(AllExp);  %for all experiments
for roi=1:LE  
Exp=strcat(init.roidirname,num2str(AllExp(roi)));
switch usr
    case 'Jacob',  expinfo=A002_JK_Condensin_with_plectonemes_expinfo(expi,AllExp(roi));
    case 'Eugene', expinfo=A002_EK_Condensin_with_plectonemes_expinfo(expi,AllExp(roi));
end
    SaveName=char(strcat(outpath0,'EKMcp_A020_',Exp)); 
    PlotSaveName=char(strcat(outpath1,'EKMcp_A020_',Exp)); 
    
        
%% Kymographs
if actions.buildkymographs     
    if mod(roi,1)==0, disp(strcat('Building kymograph:',expname,'Roi:',num2str(init.AllExp(roi)),':Exps to work through:',num2str(LE-roi+1)));end     
    %% 1)make two kymographs 
        Channel_list=[{'DNA\'}, {expinfo.labelname}];  
        dna_pth=char(strcat(generaldatapth, Exp,'\', Channel_list(1)));
        condensin_pth=char(strcat(generaldatapth, Exp,'\', Channel_list(2)));
        
        if expi~=-1
            kymo_DNA=kym_build_kymo_from_movie(dna_pth,expinfo);       
            kymo_Cnd=kym_build_kymo_from_movie(condensin_pth,expinfo); 
        else
            [kymo_DNA,kymo_Cnd]=kym_simulate_extrusion(0);
        end
        save(strcat(SaveName, '_allresults.mat'), 'kymo_DNA','kymo_Cnd');
        
        if actions.backsaving %optional back-saving a kymograph  for later use (for example, using 'boxtrack')
          dlmwrite([generaldatapth, Exp,'\kymograph\EKMcp_A020_Kymograph_MukBEF.txt'],kymo_Cnd);
        end
        
        if 1
        figure(64); 
            subplot(2,1,1); pcolor(kymo_DNA'); shading flat; colormap hot; 
            title('DNA');
            subplot(2,1,2); pcolor(kymo_Cnd'); shading flat; colormap hot; 
            title('Condensin');
            %[~]=ginput(1); 
            saveas(gcf,strcat(PlotSaveName, '_kymographs.jpg'),'jpg'); 
            pause(2);
            close(gcf);
        end
    end   

%% Analysis   
if actions.peakdetection        
      if mod(roi,1)==0, disp(strcat('Analyzing:',expname,'Roi:',num2str(roi),':Exps to work through:',num2str(LE-roi+1)));end     
  
    load(strcat(SaveName, '_allresults.mat'));
    %% 2 peak analysis on each profile of the kymographs
        kymo_DNA=kym_smooth_time_or_place(kymo_DNA,init.t_smooth,'time');
        kymo_Cnd=kym_smooth_time_or_place(kymo_Cnd,init.t_smooth,'time');  
        kymo_DNA=kym_smooth_time_or_place(kymo_DNA,init.x_smooth,'place');
        kymo_Cnd=kym_smooth_time_or_place(kymo_Cnd,init.x_smooth,'place');  
    
    
    %get some general properties  
        %% Condensin precook
        levels_Cnd=kym_get_signal_levels(kymo_Cnd,init);
        kymo_Cnd_peaks=kym_keep_only_peaks(kymo_Cnd,levels_Cnd); %shave off loops  
        tresH_Cnd=6*levels_Cnd.noise_all;
     
       
        %% DNA precook
        %get a first avalaution of raw DNA signal (background etc.)
        [levels_DNA,~,~]=kym_get_signal_levels_DNA(kymo_DNA,init,0);  
        %convert_DNA counts to_genomic_percentage:
        kymo_DNA=kym_convert_dna_kymo(kymo_DNA,levels_DNA);  
        %repeat levels  now with percentages:
        [levels_DNA,kymo_DNA_hat,kymo_DNA_peaks]=kym_get_signal_levels_DNA(kymo_DNA,init,0);  
        tresH_DNA=6*levels_DNA.level.noise_all;
          
        %% Positions DNA and Condensin.
        %Note that the residu is the DNA tether for both, in % of tether
        %total       
        info_DNA=kym_peakfitperkymographline(kymo_DNA_peaks,kymo_DNA_hat, 'peeling_and_clustering', psf_est,tresH_DNA);   %build spot info!        
        info_Cnd=kym_peakfitperkymographline(kymo_Cnd_peaks,kymo_DNA_hat,'just_treshold', psf_est,tresH_Cnd); %build spot info!
            
        save(strcat(SaveName,'_allresults.mat'),'info_DNA','info_Cnd',...
                            'kymo_DNA_peaks', 'kymo_Cnd_peaks', 'levels_DNA','levels_Cnd', '-append');    
end

if actions.smallpostprocessing
     Exp=strcat('ROI',num2str(AllExp(roi)));  
     disp(strcat('Building&saving plots: Exps to work through:',num2str(LE-roi)));
     load(strcat(SaveName, '_allresults.mat')); 
     [info_DNA,info_Cnd]=kym_get_length_and_densities(info_DNA,info_Cnd,kymo_DNA); 
      save(strcat(SaveName,'_allresults.mat'),'info_DNA','info_Cnd', '-append');
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
        
           figure(145); 
            subplot(3,1,1); pcolor(kymo_DNA'); shading flat; colormap hot; hold on;
            plot( info_DNA.pos_frameno,info_DNA.pos_X_subpix+expinfo.channelshift, 'bo','Markersize',3); hold on;           
            ylim([0 length(kymo_DNA(1,:))]);
            xlim([0 length(kymo_DNA(:,1))]);
            title('DNA');
            subplot(3,1,2); pcolor(kymo_Cnd'); shading flat; colormap hot; hold on;
            plot(info_Cnd.pos_frameno,info_Cnd.pos_X_subpix,  'bo','Markersize',3);
            ylim([0 length(kymo_DNA(1,:))]);
            xlim([0 length(kymo_DNA(:,1))]);
            title('Condensin');
            pause(0.5);
            %[~]=ginput(1); 
            subplot(3,1,3); 
            plot(info_DNA.pos_frameno,info_DNA.pos_X_subpix+expinfo.channelshift,  'bo','Markersize',4); hold on;
            plot(info_Cnd.pos_frameno, info_Cnd.pos_X_subpix, 'ro','Markersize',2);
            ylim([0 length(kymo_DNA(1,:))]);
            xlim([0 length(kymo_DNA(:,1))]);
            title(Exp);
            %legend('DNA','Condensin')
            %[sh,~,~]=ginput(2); 
            %xshift=sh(2)-sh(1) %output shift
            
            saveas(gcf,strcat(PlotSaveName, '_kymographs_overlays.jpg'),'jpg'); 
            pause(0.5); close(gcf);
            figure(169); 
            subplot(2,1,1); pcolor(kymo_DNA'); shading flat; colormap hot; hold on;
            plot( info_DNA.pos_frameno,info_DNA.pos_X_subpix+expinfo.channelshift, 'bo','Markersize',3); hold on;           
            ylim([0 length(kymo_DNA(1,:))]);
            xlim([0 length(kymo_DNA(:,1))]);
            subplot(2,1,2);
            plot(info_DNA.pos_frameno,info_DNA.content_clustercont, 'ro','Markersize',2); hold on;           
            plot(sum(kymo_DNA_peaks'),'b');
            plot(sum(kymo_DNA_hat'), 'g-');
            xlabel('frame');
            ylabel('intensity, %');
            %legend('local per peak','peak profile sum','tether sum');
            saveas(gcf,strcat(PlotSaveName, '_content_of_peaks.jpg'),'jpg');  
            pause(0.5); close(gcf);
            dum=1;
            
end
end


  



        