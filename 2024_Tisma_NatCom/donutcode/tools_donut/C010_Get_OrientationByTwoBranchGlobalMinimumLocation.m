 function Aligned=C010_Get_OrientationByTwoBranchGlobalMinimumLocation(Chromosome,Aligned,initval);
     %This function uses the density curves of two pre-detrmined bracnhes
     %to infer the orientation of the cell
     
    Aligned.Orientation='Unknown';  
    Eval_A=nanmin(Chromosome.PolarContourContent(Aligned.Orig.StartStopIndices));  %minimum first branch
    Eval_B=nanmin(Chromosome.PolarContourContent(Aligned.Orig.StopStartIndices)); %minimum second branch
        
    %determine the orientation, depending on strain type  
    if ((Eval_A<=Eval_B)&strcmp(initval.straintype,'type 1')) %Correct type 1: branch A contains gap
        HeadsOrTails='Heads';
    end   
    if ((Eval_A>Eval_B)&strcmp(initval.straintype,'type 1')) %Flipped type 1: branch B contains gap
        HeadsOrTails='Tails';
    end    
    if ((Eval_A>Eval_B)&strcmp(initval.straintype,'type 2')) %Correct type 2: branch B contains gap
        HeadsOrTails='Heads';
    end
    if ((Eval_A<=Eval_B)&strcmp(initval.straintype,'type 2')) %Flipped type 2: branch A contains gap
        HeadsOrTails='Tails';
    end
    
    Aligned.Orientation=HeadsOrTails;
    %Flip the indices if necessary      
    if (strcmp(Aligned.Orientation,'Tails')&&initval.FlipOriention);  %i ncorrect orientation; flip it!
        buf=Aligned.Orig.StartStopIndices;
        Aligned.Orig.StartStopIndices=fliplr(Aligned.Orig.StopStartIndices);  %swap&flip
        Aligned.Orig.StopStartIndices=fliplr(buf);                       %swap&flip 
        Aligned.Orig.AllIndices=[Aligned.Orig.StartStopIndices Aligned.Orig.StopStartIndices(2:end-1)];
      
    end
    
    