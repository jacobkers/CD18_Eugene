function [Aligned,ResampledPropertyCurveOrMatrix]=C025_FromEquiAngle2EquiDist(Aligned,CurveOrMatrix,CellName);
    %This function converts an original chromosome property vs. equidistant angle in
    %curves of fixed length with equidistant point
    %By definition, output total distance=100; sum of content=100;
    %JacobKers 2016 
    Ipol_Axis=linspace(0,100,100);  %by definition
    [rr,cc]=size(CurveOrMatrix);
        
    %% 3) CLEANING
    %Sort the indices (necessary for interpolation ec.). Result is a
    %continuous axis 0....1...2 where 1 is at the StopLabel position, 0 and 2
    %depict StartLabel.     
    
    NormDistAxis=Aligned.Orig.NormDist;    
    sel=find(~isnan(NormDistAxis));             %check for NaNs
    
    CurveOrMatrix_n=[];
    for ii=1:rr, CurveOrMatrix_n(ii,:)=CurveOrMatrix(ii,sel); end
    NormDistAxis=NormDistAxis(sel);          
  
    CurveOrMatrix_s=[];
    [NormDistAxis,idx]=sort(NormDistAxis);  %to be sure
    for ii=1:rr, CurveOrMatrix_s(ii,:)=CurveOrMatrix_n(ii,idx); end
    
    CurveOrMatrix_u=[];
    [NormDistAxis,idx2]=unique(NormDistAxis);  %remove doubles
    for ii=1:rr, CurveOrMatrix_u(ii,:)=CurveOrMatrix_s(ii,idx2); end
 
    CurveOrMatrix=CurveOrMatrix_u;
    
    NormDistAxis=100*(NormDistAxis-min(NormDistAxis))...
                  /(max(NormDistAxis)-min(NormDistAxis)); %re-normalize to be sure
   
    %% 5 EQUIDISTANCE : Up to here, all analysis is still versus the original annular axis.
    % Now, Interpolate the contents vs. outer contour length axis
    Ipol_CurveOrMatrix=0*repmat(Ipol_Axis,rr,1);
    for ii=1:rr
        if length(NormDistAxis)>2        
            Ipol_CurveOrMatrix(ii,:)=interp1(NormDistAxis,CurveOrMatrix(ii,:),Ipol_Axis);
        else
            Ipol_CurveOrMatrix(ii,:)=NaN*Ipol_Axis;
        end  
    end    
    
    %%  NORMALIZATION  
    % 6 normalize content of branches on whole Chromosome. Note that this is still
    %integrated vs. annular sections        
    Aligned.Dist.NormAxis=Ipol_Axis;
    ResampledPropertyCurveOrMatrix=Ipol_CurveOrMatrix;
    
    
       