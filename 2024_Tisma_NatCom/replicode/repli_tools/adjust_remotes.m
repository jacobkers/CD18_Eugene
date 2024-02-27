function outpth=adjust_remotes(inpth)     
    
    mapped_1a='N:\tnw\BN\CD\Shared\'; Lm1a=length(mapped_1a);
    mapped_1b='X:\tnw\BN\Shared\'; Lm1b=length(mapped_1b);    
    mapped_2a='M:\tnw\bn\cd\Shared\'; Lm2a=length(mapped_2a);
    mapped_2b='V:\tnw\bn\cd\Shared\'; Lm2b=length(mapped_2b);   
    outpth=inpth;    
    if contains(inpth, mapped_1a) && exist(mapped_1b,'dir')==7 
          outpth=[mapped_1b inpth(Lm1a+1:end)];
    end
    if contains(inpth, mapped_1b) && exist(mapped_1a,'dir')==7 
          outpth=[mapped_1a inpth(Lm1b+1:end)];
    end    
    if contains(inpth, mapped_2a) && exist(mapped_2b,'dir')==7 
          outpth=[mapped_2b inpth(Lm2a+1:end)];
    end
    if contains(inpth, mapped_2b) && exist(mapped_2a,'dir')==7 
          outpth=[mapped_2a inpth(Lm2b+1:end)];
    end