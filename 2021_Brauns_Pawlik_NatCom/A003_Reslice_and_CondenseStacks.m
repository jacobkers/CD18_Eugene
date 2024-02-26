function A003_Reslice_and_CondenseStacks
%JWJK_A:-------------------------------------------------------------------
%Title: cleaning of images
%
%Summary: this function tilts patterns so that each frame represents a X-T kymograph
%Next, these are condensed in XY squares so the effect is, that one has a
%collection of time-boxes. These can be used for top-bottom correlation
%work (such as done with A510_Top_and_Bottom_Correlation_XT)
%
%Input: 
%Output:

%Reference: Project: BN_CD16_Greg ; Jacob Kers 2016; Github 
%:JWJK_A-------------------------------------------------------------------

close all;
tiles=70; %tiles x tiles image
%% some setting up
close all;
if ismac, DirSep='/';else DirSep='\';end;
addpath(genpath(pwd)); 
  %specifying zoomin on the fourier plot
switch 3   
  case 2       
        initval.mainpath='D:\jkerssemakers\_Data\CD\2016_Greg\20180903 Correlation\movies_test_cln\';
        initval.searchtemplate='*.tif';
  case 3       
        initval.mainpath='D:\jkerssemakers\_Data\CD\2016_Greg\20181009_bottomtop_cleaned\';
        initval.searchtemplate='*.tif';
  end
  
    
%% build a result directory (or overwrite it)
    ReslicedSaveDir=strcat(initval.mainpath,'Resliced',DirSep);

    if ~isdir(ReslicedSaveDir)
        mkdir(ReslicedSaveDir);  
    end
    

    
%% find all files in the subdirectory that contain the 'filestring' template
    %and return a list of them (including names, path)
    impth=strcat(initval.mainpath);
    stacklist=Scroll_ImageDirs(impth,initval.searchtemplate);
    [~,ST]=size(stacklist);
    
%% for each of these files:    
    for st_ii=1:ST   
        disp(strcat(num2str(ST-st_ii+1), 'stacks to go'));
        dirname=stacklist(st_ii).dirname;
        filname=stacklist(st_ii).filname        
        info = imfinfo(strcat(dirname,DirSep,filname));    
        [ff,~]=size(info); 
        
        
    %% step 1: load first image to get an idea of the size
        firstim=double(imread(strcat(dirname,DirSep,filname),'Index',1));
        [rr,cc]=size(firstim);
        padval=median(firstim(:));

     %% step 2: load all images, build a stack.
     
        disp('Getting background image...');
        tic
        movie=zeros(tiles,tiles,ff);
         for ii=1:ff
            %disp(ff-ii)
            load_im=double(imread(strcat(dirname,DirSep,filname),'Index',ii));;
            [BackIm,TileVals]=B005_GetBackgroundSquares(load_im,tiles);
            movie(:,:,ii)=TileVals;
         end
         
                    
         
    %% step 3: save tilted stacks for x and y
    Tilt_and_Save(movie,filname,ReslicedSaveDir,'x');
    Tilt_and_Save(movie,filname,ReslicedSaveDir,'y');

    end
    

    function Tilt_and_Save(movie,filname,ReslicedSaveDir,tiltdirection)
        
    if ismac, DirSep='/';else DirSep='\';end;
    OutFileName=strcat(filname(1:end-4),'_resliced_',tiltdirection,'.tif');
    OutPathFileName=strcat(ReslicedSaveDir,DirSep,OutFileName);
    switch tiltdirection
        case 'x',     reslicemovie=shiftdim(movie,1);
        case 'y',     reslicemovie=shiftdim(movie,2);
    end
    disp('saving resliced stack...');
    [~,~,dd]=size(reslicemovie);
    for ii=1:dd
        im=reslicemovie(:,:,ii); 
        %write
        imout=uint16(im-1);         
        if ii==1        
            imwrite(imout,OutPathFileName,'tif'); 
        else
            imwrite(imout,OutPathFileName,'WriteMode','append');
        end                    
    end
