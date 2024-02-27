function A060_Colocalization(initval)
%JWJK_A:-------------------------------------------------------------------
%Description: compare channels pair-wise and determine degree of overlap of signal 

%input: .mat database from A010/013 analysis.

%output: .mat database labeled A060

%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------


close all;
showandwait=0;
save_jpg=1;

if nargin<1
    usr='Jacob', batchrunindex=0;
    initval=A000_Repli_Init(batchrunindex,usr);
end

imoutdir=strcat(initval.pth_repli,'A060_tresholds',initval.DirSep);  
if isdir(imoutdir)
    rmdir(imoutdir,'s'); 
end
mkdir(imoutdir);   

Psf_meas=initval.Psf_est;
allframes=length(initval.Cell_Labels);

%% Initialize databases here
Allcells_Mukbef_Chro_data=struct('area_A',[]);
Allcells_Mukbef_Ori_data=struct('area_A',[]);
Allcells_Mukbef_Ter_data=struct('area_A',[]);
Allcells_Ori_Ter_data=struct('area_A',[]);
Allcells_Ori_Chro_data=struct('area_A',[]);

Allcells_Mukbef_Chro_images=struct('dummy',[]);
Allcells_Ori_Chro_images=struct('dummy',[]);
Allcells_Ori_Ter_images=struct('dummy',[]);

%% run the cells
UsedCellNames=[]; ExcludedCellNames=[]; used_jj=0;
for jj=1:allframes 
    cellno=char(initval.Cell_Labels{jj}); 
    CellName=strcat('ResultsOfCell',cellno); 
    notexcluded=CheckUserExclusion(cellno, initval);
    if notexcluded
    disp(strcat('Program:A60_experiment:',initval.expi,':',CellName,'ColocalizationAnalysis..', num2str(allframes-jj+1), 'cells to go'));  
    used_jj=used_jj+1;
    pic=[];
    %first, splits image into spots
    
    for ii=1:initval.numberofchannels
        [refpic,cellmask]=Get_channel_pic(initval,cellno,2); %chromosome
        [pic,~]=Get_channel_pic(initval,cellno,ii);
        if ii==1
            [rr,cc]=size(pic);
            picstack=zeros(rr,cc,initval.numberofchannels);
        end
        picstack(:,:,ii)=pic;
        if showandwait
            figure(1);
            subplot(2,3,ii+1);pcolor(pic); shading flat, colormap hot; axis equal; axis tight;
        end
    end
    
    if showandwait
        figure(1);
        figure(1); subplot(2,3,1);pcolor(cellmask); shading flat, colormap hot; axis equal; axis tight;    
    end
    
    [~,~,chan_no]=size(picstack);
    %set the channels
    Chro=picstack(:,:,2)';      
    Ori=picstack(:,:,3)';
    Ter=picstack(:,:,4)';
    %set the tresholds
    if isfield(initval, 'tresholdperchannel')
         tr_dna=initval.tresholdperchannel(2);
         tr_ori=initval.tresholdperchannel(3);
         tr_ter=initval.tresholdperchannel(4);
         tr_mukB=initval.tresholdperchannel(5);  
    else       %auto-settings 
        tr_dna=NaN;
        tr_ori=NaN;
        tr_ter=NaN;
        tr_mukB=NaN;
    end
    [Allcells_Ori_Chro_data,Allcells_Ori_Chro_images]=Localize_2Channels(Ori,Chro,Allcells_Ori_Chro_data,[tr_ori tr_dna],used_jj,showandwait,Allcells_Ori_Chro_images);
    [Allcells_Ori_Ter_data,Allcells_Ori_Ter_images]  =Localize_2Channels(Ori,Ter,Allcells_Ori_Ter_data,[tr_ori tr_ter],used_jj,showandwait,Allcells_Ori_Ter_images);
    
    if chan_no>4
        MukB=picstack(:,:,5)';  
        [Allcells_Mukbef_Chro_data,Allcells_Mukbef_Chro_images]=Localize_2Channels(MukB,Chro,Allcells_Mukbef_Chro_data,[tr_mukB tr_dna],used_jj,0,Allcells_Mukbef_Chro_images);
        [Allcells_Mukbef_Ori_data,~]=Localize_2Channels(MukB,Ori,Allcells_Mukbef_Ori_data,[tr_mukB tr_ori],used_jj,showandwait);
        [Allcells_Mukbef_Ter_data,~]=Localize_2Channels(MukB,Ter,Allcells_Mukbef_Ter_data,[tr_mukB tr_ter],used_jj,showandwait);
    end
    
    UsedCellNames=[UsedCellNames;
                   {cellno}];
               
               
    %% build a plot: collect data we need
    %note some swapping coordinates etc. to get it right
    im_mukb=picstack(:,:,5)';
    im_dna=picstack(:,:,2)';
    im_ori=picstack(:,:,3)';
    im_ter=picstack(:,:,4)';
    
    im_tres_mukb=Allcells_Mukbef_Chro_images(used_jj).masked_image_A;
    im_tres_dna=Allcells_Mukbef_Chro_images(used_jj).masked_image_B;
    im_tres_ori=Allcells_Ori_Chro_images(used_jj).masked_image_A;
    im_tres_ter=Allcells_Ori_Ter_images(used_jj).masked_image_B;
     
     set(figure(1), 'visible','off');
        %check xy-pix etc
        subplot(2,4,1); pcolor(im_dna); shading flat; colormap jet; %dna 
        title('DNA');
        axis equal;
        subplot(2,4,2); pcolor(im_ori); shading flat; colormap jet; %ori
        title('Ori'); axis equal;
        subplot(2,4,3); pcolor(im_ter); shading flat; colormap jet; %ter
        title('Ter'); axis equal;
        subplot(2,4,4); pcolor(im_mukb); shading flat; colormap jet; %mukb
        title('MukB'); axis equal;      
        subplot(2,4,5);
         pcolor(im_tres_dna); shading flat; colormap jet;        
         axis equal;
         subplot(2,4,6);
         pcolor(im_tres_ori); shading flat; colormap jet;
         axis equal;
         subplot(2,4,7);
         pcolor(im_tres_ter); shading flat; colormap jet;
         axis equal;
         
       subplot(2,4,8);
         pcolor(im_tres_mukb); shading flat; colormap jet;
         axis equal;
                   
         outname=[CellName, '_A060_tresholds_overview']; 
         
         if save_jpg
            saveas(gcf,[imoutdir, outname, '_channels.jpg' ]);
            
         end
        close all;
    
               
               
    
    else
    disp(['excluded:' CellName]);
    ExcludedCellNames=[ExcludedCellNames;
                   {cellno}];
    end
end

save([initval.pth_repli, '\',initval.expi, '_A060_AllCellsResults.mat'], ...
                     'Allcells_Ori_Ter_data');            
if chan_no>4
    save([initval.pth_repli, '\',initval.expi, '_A060_AllCellsResults.mat'], ...
                     'Allcells_Mukbef_Chro_images',...
                     'Allcells_Ori_Chro_images',...
                     'Allcells_Ori_Ter_images',...
                     'Allcells_Mukbef_Chro_data',...
                     'Allcells_Mukbef_Ori_data',...
                     'Allcells_Mukbef_Ter_data', 'UsedCellNames', 'ExcludedCellNames', '-append');            
end
                 
                 
 function [Allcells_ChanA_ChanB_data,Allcells_ChanA_ChanB_images]=Localize_2Channels(ChanA,ChanB,Allcells_ChanA_ChanB_data,ABtresholds,jj,showandwait,Allcells_ChanA_ChanB_images);
    %compare two channels and store the results
    analysistype='bytresholding';
      
    switch analysistype 
        case 'bytresholding'  
            results=Localizer_bytresholding(ChanA,ChanB, showandwait,ABtresholds);
            %collect some quantities
            Allcells_ChanA_ChanB_data.area_A(jj,1)=results.area_A;
            Allcells_ChanA_ChanB_data.intensity_A(jj,1)=results.intensity_A;
            Allcells_ChanA_ChanB_data.area_B(jj,1)=results.area_B;
            Allcells_ChanA_ChanB_data.intensity_B(jj,1)=results.intensity_B;
            Allcells_ChanA_ChanB_data.area_AB(jj,1)=results.area_AB;
            Allcells_ChanA_ChanB_data.intensityA_inAB(jj,1)=results.intensityA_inAB;
            Allcells_ChanA_ChanB_data.intensityB_inAB(jj,1)=results.intensityB_inAB;
            Allcells_ChanA_ChanB_data.area_perc_AB_over_A(jj,1)=results.area_perc_AB_over_A;
            Allcells_ChanA_ChanB_data.area_perc_AB_over_B(jj,1)=results.area_perc_AB_over_B;
            Allcells_ChanA_ChanB_data.intensity_percA_AB_over_A(jj,1)=results.intensity_percA_AB_over_A;
            Allcells_ChanA_ChanB_data.intensity_percB_AB_over_B(jj,1)=results.intensity_percB_AB_over_B;
            
            if nargin>5  %store images    
            Allcells_ChanA_ChanB_images(jj).masked_image_A=results.masked_image_A;
            Allcells_ChanA_ChanB_images(jj).masked_image_B=results.masked_image_B;;
            Allcells_ChanA_ChanB_images(jj).masked_image_A_AND_B=results.masked_image_A_AND_B;
            dum=1;
            else
                Allcells_ChanA_ChanB_images(jj).dummy=[];
            end
    end
    
