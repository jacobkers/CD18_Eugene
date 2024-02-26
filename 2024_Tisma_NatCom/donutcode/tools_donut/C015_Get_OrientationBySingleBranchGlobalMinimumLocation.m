function Aligned=C015_Get_OrientationBySingleBranchGlobalMinimumLocation(Chromosome,Aligned,initval);
     %This function uses the density curves of two pre-detrmined bracnhes
     %to infer the orientation of the cell
     
    Aligned.Orientation='Unknown';  
    demoplot=0;
        
    [Eval,idx]=nanmin(Chromosome.PolarContourContent(Aligned.Orig.AllIndices));  
    %index points to place of global minimum starting from startlabel 
    %globalMinIndex=Aligned.Orig.AllIndices(idx);
    globalMinIndex=idx;
    MarkerIndex=find(Aligned.Orig.MarkerLabel.OriIndex==Aligned.Orig.AllIndices);
    %index points to place of marker starting from startlabel 
    
    if MarkerIndex>1
        MinComesFirst=globalMinIndex<MarkerIndex; 
    else
        MinComesFirst=globalMinIndex>max(Aligned.Orig.AllIndices)/2; 
    end
    
    %determine the orientation, depending on strain type  
    if (MinComesFirst&strcmp(initval.straintype,'type 1')) %Correct type 1: branch A contains gap
        HeadsOrTails='Heads';
    end   
    if (~MinComesFirst&strcmp(initval.straintype,'type 1')) %Flipped type 1: branch B contains gap
        HeadsOrTails='Tails';
    end    
    if ((~MinComesFirst)&strcmp(initval.straintype,'type 2')) %Correct type 2: branch B contains gap
        HeadsOrTails='Heads';
    end
    if ((MinComesFirst)&strcmp(initval.straintype,'type 2')) %Flipped type 2: branch A contains gap
        HeadsOrTails='Tails';
    end
 
    Aligned.Orientation=HeadsOrTails;
    
    %plot menu
    if demoplot
        subplot(2,1,1);
        plot(Chromosome.PolarContourContent(Aligned.Orig.AllIndices),'k-'); hold on;
        stem(MarkerIndex,Chromosome.PolarContourContent(MarkerIndex),'bo');
        stem(globalMinIndex,Chromosome.PolarContourContent(globalMinIndex),'ko');
        text(globalMinIndex+10,Chromosome.PolarContourContent(globalMinIndex),HeadsOrTails);
        title('before correction');
        legend('contour','marker','minimum');
        hold off;
    end
       
    %Flip the indices if necessary      
    if (strcmp(Aligned.Orientation,'Tails')&&initval.FlipOriention);  %incorrect orientation; flip it!                     %swap&flip 
        Aligned.Orig.AllIndices=fliplr(Aligned.Orig.AllIndices);
        Aligned.Orig.StartStopIndices=Aligned.Orig.AllIndices;
        Aligned.Orig.StopStartIndices=[];  
    end
    
    if demoplot %check correction
        [Eval,idx]=nanmin(Chromosome.PolarContourContent(Aligned.Orig.AllIndices));  
        %index points to place of global minimum starting from startlabel 
        globalMinIndex=idx;
        MarkerIndex=find(Aligned.Orig.MarkerLabel.OriIndex==Aligned.Orig.AllIndices);
        %index points to place of marker starting from startlabel      
        
        subplot(2,1,2);
        plot(Chromosome.PolarContourContent(Aligned.Orig.AllIndices),'k-'); hold on;
        stem(MarkerIndex,Chromosome.PolarContourContent(MarkerIndex),'bo');
        stem(globalMinIndex,Chromosome.PolarContourContent(globalMinIndex),'ko');
        title('before corrected');
        legend('contour','marker','minimum');
        hold off;
        [~]=ginput(1);
    end
    
    
    
    
    
    
    