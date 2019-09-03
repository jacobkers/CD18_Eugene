function [pic_list,filleddir_list]=scroll_image_dirs(impth,SearchTemplate)
%JWJK_C*:-------------------------------------------------------------------
%Title: Scroll image directories
%Summary:This function generates a list of full_pathfilenames and image names 
%of sub-directories containing the right template string; to be used for
%expanded (image) data directories
%Input: source path, search template (*.[text]);
%Output:
    %1) structure list of pairs of full_pathfilenames and image names 
    %containing the right template string
    %2) list of directories that contain the specified template
%References: written by Jacob Kers, 2010 or so
%:JWJK_C*-------------------------------------------------------------------
if nargin <2
    impth=pwd;
    SearchTemplate= '*.m';
end

    pic_list=struct('dirname',[],'filname',[]);
    filleddir_list=struct('dirname',[]);
    progdir=pwd;    
    if ismac, PathSep=':';else PathSep=';';end;
    p = genpath(char(impth)); 
    %String cotaining all (sub)dirs
    lp=length(p);
    dircount=0;
    filleddircount=0;
    filcount=0;
    oldi=1; nwi=oldi;
    while dircount<lp                        %separate out directory names
        dircount=dircount+1;
        nwi=nwi+1;       
        str=char(p(dircount));
        if strcmp(str,PathSep)              %This is one directory name
            dirnm=p(oldi:nwi-2);            
            cd(dirnm);                       %read the proper file names in this directory
            FilNames=dir(SearchTemplate);
            lm=length(FilNames); 
            if lm>0
                filleddircount=filleddircount+1;
                filleddir_list(filleddircount).dirname=dirnm;
            end
            for j=1:lm 
                filcount=filcount+1;
                pic_list(filcount).dirname=dirnm;
                pic_list(filcount).filname=FilNames(j).name ;      
            end
            oldi=nwi;
            cd(progdir);
        end  
    end  
    
if nargin <2
    filleddir_list
    pic_list
end