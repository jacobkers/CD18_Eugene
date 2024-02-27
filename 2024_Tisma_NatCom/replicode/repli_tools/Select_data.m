function Cell_Labels=Select_data(initval);  
    cd(initval.maindatapath);
    FileLabels=dir(strcat(initval.searchlabel,'1*.mat'))';  %select by first channel
    CC=length(FileLabels);
    SL=length(initval.searchlabel);
    cropit=min([CC initval.shortset]);        
    for ii=1:cropit  
        Cell_Labels{ii}=FileLabels(ii) .name(SL+3:end-4);  
            %this is a unique label per cell, per frame
    end
    cd(initval.codepth);
    
    