function  [Ipol_BPAxis,Ipol_BPContent,Ipol_PropertyCurve]=C035_Get_AnyPropertyVsMaxBasePairCurve(Aligned,PropertyCurve,initval);    
     %build a BPaxis curve based on Interpolated and Normalized distances; adapt the desired property accordingly; 
     %Optionally Correct BP axis for optional two-branch site pinning

    %% 2) Build  a Basepair-axis, optionally corrected. Still in sampling units of distance
    BPAxis=cumsum(Aligned.ViaPeakLineDistance.Ipol_PeakValVsDistance);    
    if ~Aligned.OneBranch
        %What to expect:
        switch initval.straintype
         case 'type 2'
         section1=round((100-initval.StartLabelpos)+initval.StopLabelpos);  %CW1 branch length
         section2=100-section1;                         %CW2 branch length
         case 'type 1'
         section1=round(initval.StartLabelpos-initval.StopLabelpos);        %CW1 branch length
         section2=100-section1;                         %CW2 branch length
        end       
        PinIndex=Aligned.ViaPeakLineDistance.MarkerLabelIndex;  %since in such cases this is the same as the 'Stoplabel'       
        BPPin_Unadjusted=BPAxis(PinIndex);
        CorFactor=section1/BPPin_Unadjusted;
        
        if (CorFactor>1.5)|(CorFactor<0.7)
            CorFactor=1;  %do not correct
        end  
        corFactorA=linspace(1,CorFactor,PinIndex);
        corFactorB=linspace(CorFactor,1,100-PinIndex);
        BPAxis_Cor=BPAxis.*[corFactorA corFactorB];
        
       BPAxis=BPAxis_Cor;
    end
    
    Content=diff([0 BPAxis]);   %adjusted content
    
        %% 3) CLEANING
    sel=find(~isnan(BPAxis));             %check for NaNs
    
    BPAxis=BPAxis(sel);
    Content=Content(sel);
    PropertyCurve=PropertyCurve(sel); 
               
    [BPAxis,idx]=sort(BPAxis);  %to be sure
    Content=Content(idx);
    PropertyCurve=PropertyCurve(idx); 
    
    [BPAxis,idx2]=unique(BPAxis);  %remove doubles
    Content=Content(idx2);
    PropertyCurve=PropertyCurve(idx2);   
    
    BPAxis=100*(BPAxis-min(BPAxis))...
                  /(max(BPAxis)-min(BPAxis)); %re-normalize to be sure
    
    %% 5 EQUIDISTANCE : Up to here all analysis is still versus the original annular axis.
    % Now, Interpolate the contents vs. BP content    
    Ipol_BPAxis=linspace(0,100,100);  %by definition
    if length(BPAxis)>2
        Ipol_BPContent=interp1(BPAxis,Content,Ipol_BPAxis);
        Ipol_PropertyCurve=interp1(BPAxis,PropertyCurve,Ipol_BPAxis);
    else
        Ipol_BPContent=NaN*Ipol_BPAxis;
        Ipol_PropertyCurve=NaN*Ipol_BPAxis;
    end   
    if 0
            plot(Ipol_BPAxis,Ipol_BPContent); hold on;
            [~]=ginput(1);
            close(gcf);
   end
 
