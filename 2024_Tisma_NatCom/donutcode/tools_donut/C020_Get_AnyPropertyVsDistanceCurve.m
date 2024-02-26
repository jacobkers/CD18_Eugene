function [Ipol_Axis,Ipol_Curve,Ipol_NormCurve]=C020_Get_AnyPropertyVsDistanceCurve(Chromosome,Aligned,PropertyCurve,initval);
    %This function converts an original chromosome property vs. distance in
    %curves of fixed length with equidistant point
    %By definition, output total distance=100; sum of content=100;
    %JacobKers 2016

    %1)LENGTH: Get  length along the outer contour, starting out from StartLabel;      
    p1x=Chromosome.CartesianContourEdge_X(Aligned.Ori.AllIndices);
    p1y=Chromosome.CartesianContourEdge_Y(Aligned.Ori.AllIndices);
    DistAxis=[0 [cumsum(((p1x(2:end)-p1x(1:end-1)).^2+(p1y(2:end)-p1y(1:end-1)).^2).^0.5)]'];
    NormDistAxis=100*DistAxis/max(DistAxis);
    
    
    
    
    %% 3) CLEANING
    %Sort the indices (necessary for interpolation ec.). Result is a
    %continuous axis 0....1...2 where 1 is at the StopLabel position, 0 and 2
    %depict StartLabel.
     
    Curve=PropertyCurve(Aligned.Ori.AllIndices);
    
    sel=find(~isnan(NormDistAxis));             %check for NaNs
    Curve=Curve(sel);   
    NormDistAxis=NormDistAxis(sel);          
  
    [NormDistAxis,idx]=sort(NormDistAxis);  %to be sure
    Curve=Curve(idx);    
    
    [NormDistAxis,idx2]=unique(NormDistAxis);  %remove doubles
    Curve=Curve(idx2);  
    
    NormDistAxis=100*(NormDistAxis-min(NormDistAxis))...
                  /range(NormDistAxis); %re-normalize to be sure
   
    %% 5 EQUIDISTANCE : Up to here, all analysis is still versus the original annular axis.
    % Now, Interpolate the contents vs. outer contour length axis
    Ipol_Axis=linspace(0,100,100);  %by definition
    Ipol_Curve=interp1(NormDistAxis,Curve,Ipol_Axis);
    
         
    %%  NORMALIZATION  
    % 6 normalize content of branches on whole Chromosome. Note that this is still
    %integrated vs. annular sections    
    BackGroundVal=nanmin(Ipol_Curve);
    SumContentVal=nansum(Ipol_Curve-BackGroundVal);   
    Ipol_NormCurve=100*(Ipol_Curve-BackGroundVal)/SumContentVal;   
    dum=1;