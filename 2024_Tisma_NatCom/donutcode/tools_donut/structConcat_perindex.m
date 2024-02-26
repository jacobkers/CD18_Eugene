 function OutStructure=structConcat_perindex(InStructureA,InStructureB); 
       %concatenates two structures
       %Jacob Kers 2021
       [rA,cA]=size(InStructureA);
       [rB,cB]=size(InStructureB);
       ii_12=0;
       for ii_1=1:cA
           ii_12=ii_12+1;
           OutStructure(ii_12)=InStructureA(ii_1);
       end
       for ii_2=1:cB
           ii_12=ii_12+1;
           OutStructure(ii_12)=InStructureB(ii_2);
       end
        