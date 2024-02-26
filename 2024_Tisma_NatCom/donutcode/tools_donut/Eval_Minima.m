 function  minima=Eval_Minima(DensityAxis,DensityCurve,BP_axis);
     %This function avaluates mininima in a curve. We make a distinction
     %between 'srong' minimum (below 0.25 of the average curve intensity) 
     %and 'weak' minima (below 0.5)
    prf=DensityCurve;   
    treslevelWeak=1.0*nanmean(prf);
    treslevelStrong=0.5*nanmean(prf);    
    prf_left=prf(1:end-2);
    prf_mid=prf(2:end-1);
    prf_right=prf(3:end);       
    selmin=1+find(((prf_mid<prf_left)&(prf_mid<prf_right))&...
                    (prf_mid<treslevelWeak));  %strong&weak minima 
    if length(BP_axis)>1       
        mincount=nansum(selmin);
        minima.Idx=selmin;
        minima.Vals=DensityCurve(selmin);
        minima.Pos=DensityAxis(selmin);
        minima.BPpos=BP_axis(selmin);
        minima.Strength=0*minima.Idx+1;
        strongones=find(DensityCurve(selmin)<treslevelStrong);
        minima.Strength(strongones)=2;
    else
        mincount=0;
        minima.Idx=NaN;
        minima.Vals=NaN;
        minima.Pos=NaN;
        minima.BPpos=NaN;
        minima.Strength=NaN;
    end