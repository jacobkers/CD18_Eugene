function [BPAxis,DistanceAxis]=F010_Get_BasePairvsDistanceRelation(initval,Aligned)
    %This function relates BP to distance for post-processing (after
    %alignment)
    dum=1;
    if ~initval.ViaPeakVals
        BPAxis=cumsum(Aligned.ViaContourDistance.Ipol_ContentVsDistance);
        DistanceAxis=1:length(BPAxis);  
     else
        BPAxis=cumsum(Aligned.ViaPeakLineDistance.Ipol_PeakValVsDistance);
        DistanceAxis=1:length(BPAxis);  
    end