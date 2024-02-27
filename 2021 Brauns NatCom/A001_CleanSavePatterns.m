function A001_CleanSavePatterns
%JWJK_A:-------------------------------------------------------------------
%Title: cleaning of images
%
%Summary: this function applies background cleaning of patterns. 
% Background correction and artefact removal was done as follows: 
% first ,  for movies, frames were corrected for fluorescence bleaching 
% by normalizing each frame on its mean intensity value. 
% This corrected for the max. 20% intensity decay over long movies. 
% Next, two correction images were processed: 
% 1) a ‘static background’-image Imstat was made by averaging out all 
% moving (wave pattern) features of the movie stack and removing any 
% residual background level. Thus, this image only contained static 
% fluorescent features such as specks, holes and scratches. 
% 2) an ‘illumination correction’ image Imillum  was made by strongly 
% smoothening out and averaging all movie images and normalizing the result 
% to its maximum. Finally, each movie image Immovie was corrected as via the
% following image operation: Imcorrected=(Immovie- Imstat)/Imillum  .  
% This way, irregularities are suppressed and wave amplitudes on the edge 
% of each image would not be underestimated compared to the amplitudes 
% in the center of the image. 
%
%Input: Paths are set at line ~~5. Movies should be stacked tif files
%Output: cleaned tif stacks, extension 'cln'

%Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
%code designed & written Jacob Kerssemakers 2016 
%:JWJK_A-------------------------------------------------------------------

close all;

%% some setting up
close all;
if ismac, DirSep='/';else DirSep='\';end;
addpath(genpath(pwd)); 
workpath='C:\Users\jkerssemakers\CD_Data_in\Sabrina\';
  %specifying zoomin on the fourier plot
switch 1   
    case 1
        initval.maininpath=[workpath, '\20201022_Min_flow\20201022_1_without_flow\'];
        initval.mainoutpath=[workpath, '\20201022_Min_flow\20201022_1_without_flow_cl\'];
        initval.searchtemplate='*.tif';
end
  
  onlydiagnose=0;
    
%% build a result directory (or overwrite it)
    CleanSaveDir=initval.mainoutpath;

    if ~isdir(CleanSaveDir)
        mkdir(CleanSaveDir);  
    end
    

    
%% find all files in the subdirectory that contain the 'filestring' template
    %and return a list of them (including names, path)
    impth=strcat(initval.maininpath);
    stacklist=Scroll_ImageDirs(impth,initval.searchtemplate);
    [~,ST]=size(stacklist);
    
%% for each of these files:    
    for st_ii=1:ST      
        dirname=stacklist(st_ii).dirname;
        filname=stacklist(st_ii).filname        
        info = imfinfo(strcat(dirname,DirSep,filname));    
        [ff,~]=size(info); 
        OutFileName=strcat(filname(1:end-4),'_cln.tif');
        OutPathFileName=strcat(CleanSaveDir,DirSep,OutFileName);
        
    %% step 1: load first image to get an idea of the size
        firstim=double(imread(strcat(dirname,DirSep,filname),'Index',1));
        [rr,cc]=size(firstim);
        padval=median(firstim(:));
        
    %% step 1a: perform bleach correction
        disp('obtaining bleach curve correction.....');
        bleachcurve=zeros(ff,1);        
            for ii=1:ff
                im=double(imread(strcat(dirname,DirSep,filname),'Index',ii));
            bleachcurve(ii)=mean(im(:));
            end
            sel=find(bleachcurve==0); 
            bleachcurve(sel)=bleachcurve(sel-1);
            bleachcurve=bleachcurve/max(bleachcurve);
            bleachcurve=smooth(bleachcurve,6);
      
     if onlydiagnose           
        figure(6);
        plot(bleachcurve); pause(0.1); hold on;  
    end   
    correctcurve=1./bleachcurve;   
    if ~onlydiagnose  
 
         
     %% step 2: load all images, build a 'dirt_image' that contains only static features
        disp('Getting background image...');
        tic
        stacked_dirt_im=zeros(rr,cc);
        stacked_illum_im=zeros(rr,cc);
        goodcount=0;
         for ii=1:ff
            %disp(ff-ii)
            raw_im=double(imread(strcat(dirname,DirSep,filname),'Index',ii));
            
            %correct darkimages
            if sum(raw_im(:))==0 
                raw_im=raw_im+padval; 
                emptyim=1;
            else
                emptyim=0;
                goodcount=goodcount+1;
            end
            
            if ~emptyim                
            raw_im=correctcurve(ii)*raw_im;      
            medvalraw=median(raw_im);
            unspiked_im=raw_im;
            %1) replace worst worst outliers by 'clipped range'  
            [flag,cleandata]=JKD1_PRF_outlier_flag(-raw_im(:),6,0.9,'positive',0); %holes
            if ~isempty(flag) 
                unspiked_im(flag==0)=-max(cleandata); 
            end                        
            [flag,cleandata]=JKD1_PRF_outlier_flag(raw_im(:),6,0.9,'positive',0); %peaks
            if ~isempty(flag) 
                unspiked_im(flag==0)=max(cleandata); 
            end 
  
             if 0 %plot or not
                 figure
                 subplot (1,2,1);pcolor(raw_im); colormap hot; shading flat; axis square
                 title('raw');
                 subplot (1,2,2);pcolor(unspiked_im);colormap hot; shading flat; axis square
                 title('unspiked');
            pause(0.01);
            [~] =ginput(1);
            close(gcf);
            end
                    
            %) split image in smooth and object part; add both
            [ObjectImage, IllumImage,ImProps]=Split_Image(unspiked_im);
            if 0%plot or not
                 figure
                 subplot (1,2,1);pcolor(IllumImage); colormap hot; shading flat; axis square
                 title('illumination');
                 subplot (1,2,2);pcolor(ObjectImage);colormap hot; shading flat; axis square
                 title('objects');
                pause(0.01);
                [~] =ginput(1);
                close(gcf);
            end
            stacked_dirt_im=stacked_dirt_im+ObjectImage;
            stacked_illum_im=stacked_illum_im+IllumImage;
            end
         end
         
         dirt_im=stacked_dirt_im/goodcount;
         dirt_im=(dirt_im-median(dirt_im(:)));
         
         
         stacked_illum_im=GetBackgroundImage(stacked_illum_im,25);
         illum_im=stacked_illum_im/max(stacked_illum_im(:));
        toc
        if 0  %plot or not
            figure; 
                subplot (2,2,1);pcolor(dirt_im); colormap hot; shading flat; axis equal
                title('average dirt');
                subplot (2,2,2);pcolor(illum_im);colormap hot; shading flat; axis equal
                title('average illumination');
                hfrr=round(rr/2); hfcc=round(cc/2); 
                subplot (2,2,3);
                    plot(dirt_im(hfrr,:),'b-'); hold on;
                    plot(dirt_im(:,hfcc),'r-'); axis square; 
                    axis square; title('dirt cross-sections');
                subplot (2,2,4);plot(illum_im(hfrr,:),'b-'); hold on;
                    plot(illum_im(:,hfcc),'r-'); axis square; 
                    axis square; title('correction cross-sections');
                    legend('illum-mid hor', 'illum-mid ver');
                    
            pause(0.1);
            [~] =ginput(1);
        end

    %% step 3: remove background and remove outliers per image
        %figure;
        diagnostix=[];
        disp('cleaning stack...');
        tic
        for ii=1:ff
            %disp(ff-ii)
                    if mod(ii,20)==0, disp(strcat(num2str(ff-ii+1),'images and',num2str(ST-st_ii+1), 'stacks to go')); end
            im=double(imread(strcat(dirname,DirSep,filname),'Index',ii)); 
            
            %correct darkimages
            if sum(im(:))==0, 
                im=im+padval; 
                emptyim=1;
            else
                emptyim=0;
                goodcount=goodcount+1;
            end
            
            if ~emptyim
                 clean_im=((im-dirt_im)./illum_im);
                [flag,cleandata]=JKD1_PRF_outlier_flag(clean_im(:),6,0.9,'positive',0);
                if ~isempty(flag), clean_im(flag==0)=max(cleandata(:)); end
                clean_im=JKD2_IM_smoothJK(clean_im,2);
                else
                    clean_im=im;
            end
            diagnostix=[diagnostix; [mean(clean_im(:)) std(clean_im(:))] std(im(:))]; 
            
 
            %write
            imout=uint16(clean_im-1);         
            if ii==1        
                imwrite(imout,OutPathFileName,'tif'); 
            else
                imwrite(imout,OutPathFileName,'WriteMode','append');
            end  
            
            if 0  %plot or not
            figure; 
                subplot (2,2,1);pcolor(im); colormap bone; shading flat; axis square; 
                title(['raw image', num2str(ii)]);
                subplot (2,2,2);pcolor(clean_im);colormap bone; shading flat; axis square; 
                title('cleaned');
                hfrr=round(rr/2); hfcc=round(cc/2); 
                subplot (2,2,3);
                    plot(im(hfrr,:),'b-'); hold on;
                    plot(im(:,hfcc),'r-'); axis square; 
                    axis square; 
                    title('raw cross-sections');
                    ylim([min(im(:)) max(im(:))]);
                subplot (2,2,4);
                    plot(clean_im(hfrr,:),'b-'); hold on;
                    plot(clean_im(:,hfcc),'r-'); axis square; 
                    axis square; 
                    title('cleaned cross-sections');
                    ylim([min(im(:)) max(im(:))]);
                    legend('hor', 'ver');
                    
            pause(0.1);
            [~] =ginput(1);
            end
            
            
            if 0
            subplot(1,3,1);pcolor(clean_im); colormap gray; shading flat; axis equal
            subplot(1,3,2); plot(diagnostix(:,1));
            subplot(1,3,3); plot(diagnostix(:,2)); hold on;
            subplot(1,3,3); plot(diagnostix(:,3),'r'); hold off;
            pause(0.011);
            end
        end
        end
        toc
    end

