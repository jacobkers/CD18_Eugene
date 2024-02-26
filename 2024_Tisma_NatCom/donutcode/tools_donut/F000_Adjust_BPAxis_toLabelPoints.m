function BPAxis=F000_Adjust_BPAxis_toLabelPoints(BPAxis,Aligned,initval);
        switch initval.straintype
         case 'type 2'
         section1=round((100-initval.StartLabelpos)+initval.StopLabelpos);  %CW1 branch length
         section2=100-section1;                         %CW2 branch length
         case 'type 1'
         section1=round(initval.StartLabelpos-initval.StopLabelpos);        %CW1 branch length
         section2=100-section1;                         %CW2 branch length
        end   

        if initval.ViaPeakVals
            PinIndex=Aligned.ViaPeakLineDistance.MarkerLabelIndex;  %since in such cases this is the same as the 'Stoplabel'       
        else
            PinIndex=Aligned.ViaContourDistance.MarkerLabelIndex;  %since in such cases this is the same as the 'Stoplabel'             
        end
        
        BPPin_Unadjusted=BPAxis(PinIndex);
        CorFactor=section1/BPPin_Unadjusted;
        
        if (CorFactor>2)|(CorFactor<0.5)
            CorFactor=1;  %do not correct
        end  
        corFactorA=linspace(1,CorFactor,PinIndex);
        corFactorB=linspace(CorFactor,1,100-PinIndex);
        BPAxis_Cor=BPAxis.*[corFactorA corFactorB];
        
        BPAxis=BPAxis_Cor;

    
     