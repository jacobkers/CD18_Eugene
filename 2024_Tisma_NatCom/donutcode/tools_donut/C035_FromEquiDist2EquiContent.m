function  [Aligned,ResampledPropertyCurveOrMatrix]=C035_FromEquiDist2EquiContent(Aligned,PropertyCurveOrMatrix,initval);    
     %build a BPaxis curve based on Interpolated and Normalized distances; adapt the desired property accordingly; 
     %Optionally Correct BP axis for optional two-branch site pinning
     %before using it as an interpolation axis
    
    Ipol_BPAxis=linspace(0,100,100);  %by definition
    [rr,cc]=size(PropertyCurveOrMatrix);
    
    
    
    
    %% 2) Build  a Basepair-axis, optionally corrected. Still in sampling units of distance
    BPAxis=Aligned.Dist.NormCumDensity; 
    MarkerPerc=Aligned.Orig.MarkerLabel.MarkerBPPerc;
    MarkerexpectedPerc=Aligned.Orig.MarkerLabel.MarkerBPPercExpected;
    [~,PinIndex]=min(abs(BPAxis-MarkerPerc));  %where we find the marker now
    
    CorrectionOK=1;
    %optional two-branch pinning
    if Aligned.OneBranch==0
        CorFactor=MarkerexpectedPerc/MarkerPerc;      
         if ((CorFactor>1.5)|(CorFactor<0.6))&(~initval.PassAllCells);
             CorFactor=NaN;  %do not correct
             CorrectionOK=0;
         end  
        %section-wise linear correction such that the markerindex
        corFactorA=linspace(1,CorFactor,PinIndex);
        corFactorB=linspace(CorFactor,1,100-PinIndex);
        CorBPperDist=[corFactorA corFactorB];
        BPAxis_Cor=BPAxis.*CorBPperDist;        
        BPAxis=BPAxis_Cor;
    end
    UncleanedBPdistAxis=BPAxis;
    if CorrectionOK         
    %% 3) CLEANING
    sel=find(~isnan(BPAxis));             %check for NaNs    
    BPAxis=BPAxis(sel);
    PropertyCurveOrMatrix_c=[];
    for ii=1:rr, PropertyCurveOrMatrix_c(ii,:)=PropertyCurveOrMatrix(ii,sel); end             
    
    PropertyCurveOrMatrix_s=[];
    [BPAxis,idx]=sort(BPAxis);  %to be sure    
    for ii=1:rr, PropertyCurveOrMatrix_s(ii,:)=PropertyCurveOrMatrix_c(ii,idx); end    
    
    PropertyCurveOrMatrix_dd=[];
    [BPAxis,idx2]=unique(BPAxis);  %remove doubles
    
    for ii=1:rr, PropertyCurveOrMatrix_dd(ii,:)=PropertyCurveOrMatrix(ii,idx2); end    
     
    PropertyCurveOrMatrix=PropertyCurveOrMatrix_dd;
    
    BPAxis=100*(BPAxis-min(BPAxis))...
                  /(max(BPAxis)-min(BPAxis)); %re-normalize to be sure
              
    
    %% 5 EQUIDISTANCE : Up to here all analysis is still versus the original annular axis.
    % Now, Interpolate the contents vs. BP content  
    Ipol_PropertyCurveOrMatrix=0*repmat(Ipol_BPAxis,rr,1);
    for ii=1:rr
        if length(BPAxis)>2
            Ipol_PropertyCurveOrMatrix(ii,:)=interp1(BPAxis,PropertyCurveOrMatrix(ii,:),Ipol_BPAxis);
        else
            Ipol_PropertyCurveOrMatrix(ii,:)=NaN*Ipol_BPAxis;
        end  
    end
        Aligned.BP.NormAxis=Ipol_BPAxis;
        %Aligned.Dist.NormCumDensityMarkercorrected=BPAxis;
        %this is the BP axis vs. equidistance, but now corrected with the
        %marker position. (If single branch, it is identical)
        Aligned.Dist.NormCumDensityMarkercorrected=UncleanedBPdistAxis;
        %[~,PinIndexCor]=min(abs(Ipol_BPAxis-MarkerPerc));  %where we find the marker now
        Aligned.BP.MarkerLabel.MarkerBPPerc=UncleanedBPdistAxis(PinIndex);    
        ResampledPropertyCurveOrMatrix=Ipol_PropertyCurveOrMatrix;
    else
        Aligned.BP.NormAxis=Ipol_BPAxis;
        Aligned.Dist.NormCumDensityMarkercorrected=NaN.*Aligned.Dist.NormAxis;
        Aligned.BP.MarkerLabel.MarkerBPPerc=NaN; 
        ResampledPropertyCurveOrMatrix=NaN*repmat(Ipol_BPAxis,rr,1);
    end
    Aligned.BP.CorrectionOK=CorrectionOK;
    
 
