function Cell_Labels=Select_data(initval);  
    cd(initval.maindatapath);
    chan1_string=[initval.searchlabel,'1'];
    FileLabels=dir(['*',chan1_string, '*.mat'])';  %select by first channel
    CC=length(FileLabels);
    SL=length(initval.searchlabel);
    cropit=min([CC initval.shortset]);        
    for ii=1:cropit  
        thislabel=FileLabels(ii).name;
        stri=strfind(thislabel,chan1_string)+SL+2;
        Cell_Labels{ii}=thislabel(stri:end-4);  
            %this is a unique label per cell, per frame
    end
    cd(initval.codepth);
    
    