function dir_list=find_data_dir(start_pth,SearchTemplate,first_one)
%find a data directory
if nargin <3
    first_one=1;
    start_pth='D:\jkerssemakers\Dropbox\CD_Data_out\';      %start path
    SearchTemplate= '20201022_Min_flow\20201022_1_without_flow\';         %directory to find
end
    dir_list=[];
    progdir=pwd;    
    PathSep=';';
    pp = genpath(char(start_pth)); 
    %String containing all (sub)dirs
    dirseps=strfind(pp,PathSep);
    Ldirs=length(dirseps);
    dirseps=[dirseps Ldirs+1];
    lp=length(pp);
    dircount=0;
    oldi=1; nwi=oldi;
    for ii=1:Ldirs                        
        start_i=dirseps(ii)+1;
        stop_i=dirseps(ii+1)-1;
        dirnm=[pp(start_i:stop_i) '\'];
        %bottomdir=dirnm(thisdirseps(end)+1:end);
        if contains(dirnm, SearchTemplate)  %contains dir
            start_i=strfind(dirnm, SearchTemplate);
            stop_i=start_i-1+length(SearchTemplate);
            if stop_i==length(dirnm) %bottom part found!
                dircount=dircount+1;
                disp(dirnm);
                dir_list{dircount}=dirnm;
            end
        end
    end  
    if ~isempty(dir_list)&& first_one==1
        dir_list=char(dir_list{1});
    end
        
  