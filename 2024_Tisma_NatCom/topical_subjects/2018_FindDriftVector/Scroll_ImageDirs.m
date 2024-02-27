function pic_list=Scroll_ImageDirs(initval)
     %This function generates a list of full_pathfilenames containing the right
     %string
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    progdir=pwd;
    impth=strcat(initval.mainpath,initval.imagesubdirpath);
    p = genpath(char(impth)); 
    %String cotaining all (sub)dirs
    lp=length(p);
    dircount=0;
    filcount=0;
    oldi=1; nwi=oldi;
    while dircount<lp                        %separate out directory names
        dircount=dircount+1;
        nwi=nwi+1;       
        str=char(p(dircount));
        if strcmp(str,';')              %This is one directory name
            dirnm=p(oldi:nwi-2);
            cd(dirnm);                       %read the proper file names in this directory
            FilNames=dir(initval.filestring);
            lm=length(FilNames);          
            for j=1:lm 
                filcount=filcount+1;
                pic_list(filcount).dirname=dirnm;
                pic_list(filcount).filname=FilNames(j).name ;      
            end
            oldi=nwi;
            cd(progdir);
        end  
    end  