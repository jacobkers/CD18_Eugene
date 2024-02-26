 function OutStructure=structConcat(InStructureA,InStructureB); 
       %concatenates all fields of two structures with identical column fields
       %Jacob Kers 2021
       emptyA=isempty(InStructureA);
       emptyB=isempty(InStructureB);
       X=[];
       if ~emptyA, X=fieldnames(InStructureA); end  
       if ~emptyB&emptyA, X=fieldnames(InStructureB); end
       
       if ~isempty(X)
          for ii=1:length(X)
             if ~emptyA, fcA=InStructureA.(X{ii});else fcA=[];end
             if ~emptyB, fcB=InStructureB.(X{ii});else fcB=[];end
             OutStructure.(X{ii})=[fcA; fcB];
          end
       else
           Outstructure=[];
       end
       dum=1;