 function Aligned=B020_AlignChromosomeContent(Chromosome,StopLabel,StartLabel,MarkerLabel,initval);
% 'Use this section for a Quicksheet'
%-------------------------------------------------------------------------
    %This function returns an adapted content curve, re-interpolated to
    %align on StopLabel and StartLabel positions (opotional). The order of
    %contour points is set by the radial axis. We start at StartLabel
    %------------------------------------------------
% 'End of Quicksheet section'
%JacobKers2016

Aligned=struct('Ori',[]);


%% Check number of label positions to use: one or two 
if StartLabel.spotPosAngle==StopLabel.spotPosAngle,
    Aligned.OneBranch=1;
else
    Aligned.OneBranch=0;
end

%1) find start(StartLabel) and stop(StopLabel) points on contour curve
    [minStopLabeldist,StopLabelIndex]=min...   %annular index of nearest point
    (((Chromosome.CartesianContourEdge_X-StopLabel.spotX).^2+...
    (Chromosome.CartesianContourEdge_Y-StopLabel.spotY).^2).^0.5);        
    [minStartLabeldist,StartLabelIndex]=min...   %annular index of nearest point
    (((Chromosome.CartesianContourEdge_X-StartLabel.spotX).^2+...
    (Chromosome.CartesianContourEdge_Y-StartLabel.spotY).^2).^0.5);
    [minSMarkerLabeldist,MarkerLabelIndex]=min...   %annular index of nearest point
    (((Chromosome.CartesianContourEdge_X-MarkerLabel.spotX).^2+...
    (Chromosome.CartesianContourEdge_Y-MarkerLabel.spotY).^2).^0.5);
    
    Aligned.Ori.StartLabelIndex=StartLabelIndex;
    Aligned.Ori.StopLabelIndex=StopLabelIndex;
    Aligned.Ori.MarkerLabelIndex=MarkerLabelIndex;

    %% 2) find angular indices in the right order to define an always clockwise
    %StartLabel-StopLabel and StopLabel-StartLabel branch (the latter maybe empty) 
    AngleNos=length(Chromosome.AnnularAxis);
    AngularIndices=1:AngleNos;  %angles run CCW!
    if ~Aligned.OneBranch
        if StartLabelIndex<StopLabelIndex %zero-StartLabel-StopLabel
            Aligned.Ori.StartStopIndices=(StartLabelIndex:StopLabelIndex);
            Aligned.Ori.StopStartIndices=([StopLabelIndex:1:AngleNos 1:1:StartLabelIndex]);
            %Note that overlapping start and end points points are removed:
            Aligned.Ori.AllIndices=[Aligned.Ori.StartStopIndices Aligned.Ori.StopStartIndices(2:end-1)];
        end
        if StartLabelIndex>StopLabelIndex  %zero-StopLabel-StartLabel 
            Aligned.Ori.StartStopIndices=([StartLabelIndex:AngleNos 1:StopLabelIndex]);   
            Aligned.Ori.StopStartIndices=(StopLabelIndex:1:StartLabelIndex);
            %Note that overlapping start and end points points are removed:
            Aligned.Ori.AllIndices=[Aligned.Ori.StartStopIndices Aligned.Ori.StopStartIndices(2:end-1)];
        end
    else       %SINGLE BRANCH         
            Aligned.Ori.StartStopIndices=[[StartLabelIndex:AngleNos] [1:StopLabelIndex-1]];           
            Aligned.Ori.StopStartIndices=[];  
            Aligned.Ori.AllIndices=Aligned.Ori.StartStopIndices;
    end
    
    if initval.AddFlippMarkers
        LP=length(Aligned.Ori.AllIndices);
        MarkIdx=Aligned.Ori.AllIndices(ceil(LP/4));  %take the index at one quarter from the chosen label start position
        lo=max([MarkIdx-2 1]); hi=min([MarkIdx+2 LP]);
        Chromosome.PolarContourPeakVal(lo:hi)=max(Chromosome.PolarContourPeakVal);
    end
    
    
    
    %% 2a Get orientation and optionally flip the indices accordingly:
    %identify which branch contains gap: this is the branch with the minimum value
    %Since we know the relative genomic locations of StartLabel, gap and StopLabel per strain, this
    %sets the orientation. Note that this depends on the experiment. Indices are shuffled optionally and accordingly 
    if ~Aligned.OneBranch        
        Aligned=C010_Get_OrientationByTwoBranchGlobalMinimumLocation(Chromosome,Aligned,initval);
    else
        Aligned=C015_Get_OrientationBySingleBranchGlobalMinimumLocation(Chromosome,Aligned,initval);
    end

   %% ----------------------------------------------------------------------
    %Obtain density curve vs equidistant contour points, normalized to 100
    PropertyCurve=Chromosome.PolarContourPeakVal;
    [Ipol_DistAxis,Ipol_Curve,Ipol_NormCurve]=C020_Get_AnyPropertyVsDistanceCurve(Chromosome,Aligned,PropertyCurve,initval);
    Aligned.ViaContourDistance.Ipol_OuterContourAxis=Ipol_DistAxis;
    Aligned.ViaContourDistance.Ipol_ContentVsDistance=Ipol_NormCurve;
      
        %% Link Marker to new axes
     PropertyCurve=Chromosome.MappingIndex;
    [Ipol_DistAxis,Ipol_MappedIndices,~]=C020_Get_AnyPropertyVsDistanceCurve(Chromosome,Aligned,PropertyCurve,initval);
     Aligned.ViaContourDistance.Ipol_MappedIndices=Ipol_MappedIndices;
     [~,idx]=min(abs(Ipol_MappedIndices-Aligned.Ori.MarkerLabelIndex));
     Aligned.ViaContourDistance.MarkerLabeldistance=Ipol_DistAxis(idx);
     Aligned.ViaContourDistance.MarkerLabelIndex=idx;

    %% Build Basepairaxis-content; starting from interpolated distance-based curves
    PropertyCurve=Aligned.ViaContourDistance.Ipol_ContentVsDistance;
    [Ipol_BPAxis,Ipol_BPContent,~]=C030_Get_AnyPropertyVsBasePairCurve(Aligned,PropertyCurve,initval);    
    Aligned.ViaContourBasePair.Ipol_BPAxis=Ipol_BPAxis;
    Aligned.ViaContourBasePair.Ipol_ContentVsBasepair=Ipol_BPContent;
    

    %% Link Marker to new axes
    PropertyCurve=Aligned.ViaContourDistance.Ipol_MappedIndices;
    [Ipol_BPAxis,~,Ipol_BPMappedIndices]=C030_Get_AnyPropertyVsBasePairCurve(Aligned,PropertyCurve,initval);    
    Aligned.ViaContourBasePair.Ipol_BPMappedIndices=Ipol_BPMappedIndices;
    [~,idx2]=min(abs(Ipol_BPMappedIndices-Aligned.Ori.MarkerLabelIndex));
    Aligned.ViaContourBasePair.MarkerLabelBP=Ipol_BPAxis(idx2);
    Aligned.ViaContourBasePair.MarkerLabelIndex=idx2;
    dum=1;
         