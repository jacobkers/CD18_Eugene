function  maxima=Eval_Maxima(DensityAxis,DensityCurve,BP_axis);
    prf=DensityCurve;    
    treslevelWeak=1*nanmean(prf);
    treslevelStrong=2*nanmean(prf);    
    prf_left=prf(1:end-2);
    prf_mid=prf(2:end-1);
    prf_right=prf(3:end);            
    selmax=1+find(((prf_mid>prf_left)&(prf_mid>prf_right))&...
            (prf_mid>treslevelWeak));        
    maxcount=nansum(selmax);
    maxima.Idx=selmax+1;
    maxima.Vals=DensityCurve(selmax);
    maxima.Pos=DensityAxis(selmax)';   
    maxima.BPpos=BP_axis(selmax);
    maxima.Strength=0*maxima.Idx+1;
    strongones=find(DensityCurve(selmax)>treslevelStrong);
    maxima.Strength(strongones)=2;
    
    if 0
        close all;
        plot(DensityAxis,DensityCurve, 'b-'); hold on;
         plot(DensityAxis(selmax),DensityCurve(selmax), 'ro'); hold on;
        [~]=ginput(1);
        close(gcf);
    end