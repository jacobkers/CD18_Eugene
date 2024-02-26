function A002_ResliceStacks
%JWJK_A:-------------------------------------------------------------------
%Title: reslicing of cleaned images
%
%Summary: this function tilts patterns: for quicker analysis, xy-t stacks
%are tilted x-t and y-t, such that each saved frame is essentially a
%kymograph
%
%Input: directory with pre-cleaned movies
%Output: per movie, an xt and a yt movie.

%Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
%code designed & written Jacob Kerssemakers 2016 
%:JWJK_A-------------------------------------------------------------------

close all;

%% some setting up
close all;
if ismac, DirSep='/';else DirSep='\';end;
addpath(genpath(pwd)); 
  %specifying zoomin on the fourier plot
switch 5   
    case 5       
        initval.maininpath=[pwd, '\testdata\movies_test_cln\'];
        initval.mainoutpath=[pwd, '\testdata\movies_test_cln_rs\'];
        initval.searchtemplate='*.tif';
 end
  
    
%% build a result directory (or overwrite it)
    ReslicedSaveDir=initval.mainoutpath;

    if ~isdir(ReslicedSaveDir)
        mkdir(ReslicedSaveDir);  
    end
    

    
%% find all files in the subdirectory that contain the 'filestring' template
    %and return a list of them (including names, path)
    impth=strcat(initval.maininpath);
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
        movie=zeros(rr,cc,ff);
         for ii=1:ff
            %disp(ff-ii)
            movie(:,:,ii)=double(imread(strcat(dirname,DirSep,filname),'Index',ii));
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
