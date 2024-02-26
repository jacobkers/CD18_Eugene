 function Chromosome=B025_AlignChromosomeContent_OneBranch(Chromosome,StopLabel,StartLabel,initval);
% 'Use this section for a Quicksheet'
%-------------------------------------------------------------------------
    %This function returns an adapted content curve, re-interpolated to
    %align from a 'StartLabel' position. To keep savings format
    %compatibel with the 'B020' two-lalbel alignment case, empty fields are
    %created for the non-existing 'stoplabel-startlabel' branch
    
    % NEW fields added in this function 
    %Branches, aligned at expected label positions
    
    %     Ipol_StartLabelStopLabel_Content
    %     Ipol_StartLabelStopLabel_PercAxis
    
    %     Ipol_StopLabelStartLabel_Content %these will stay empty
    %     Ipol_StopLabelStartLabel_PercAxis %these will stay empty
    
    %     Ipol_OuterContourAxis
    %     Ipol_OuterContourContent
    %     Orientation: 'Heads'   
            %this is detrmined by the location of the
    %     global minimum
    %     StartLabelStopLabelBranch_CW 
    %     StopLabelStartLabelBranch_CW
             %These are indices in the original curves
    %------------------------------------------------
% 'End of Quicksheet section'
%JacobKers2016

    %1) find start(StartLabel) and stop(StopLabel) points on contour curve
    [minStopLabeldist,StopLabelIndex]=min...   %annular index of nearest point
    (((Chromosome.CartesianContourEdge_X-StopLabel.spotX).^2+...
    (Chromosome.CartesianContourEdge_Y-StopLabel.spotY).^2).^0.5);        
    [minStartLabeldist,StartLabelIndex]=min...   %annular index of nearest point
    (((Chromosome.CartesianContourEdge_X-StartLabel.spotX).^2+...
    (Chromosome.CartesianContourEdge_Y-StartLabel.spotY).^2).^0.5);

    
    %% 2) find angular indices in the right order to define an always clockwise
    %StartLabel-StartLabel branch 
    AngleNos=length(Chromosome.AnnularAxis);
    AngularIndices=1:AngleNos;  %angles run CCW!    
    Chromosome.StartLabelStopLabelBranch_CW=[[StopLabelIndex:AngleNos] [1: StartLabelIndex]];
    Chromosome.StopLabelStartLabelBranch_CW=[];      

   
    %% ----------------------------------------------------------------------
    %3)Get  length axes along the branch, starting out from StartLabel
    %Clockwise A
    p1x=Chromosome.CartesianContourEdge_X(Chromosome.StartLabelStopLabelBranch_CW);
    p1y=Chromosome.CartesianContourEdge_Y(Chromosome.StartLabelStopLabelBranch_CW);
    DistAxis_A=[0 [cumsum(((p1x(2:end)-p1x(1:end-1)).^2+(p1y(2:end)-p1y(1:end-1)).^2).^0.5)]'];
    NormDistAxis=DistAxis_A/max(DistAxis_A);    
 
    %% ----------------------------------------------------    
    %4 normalize content of branch on whole Chromosome. Note that this is still
    %integrated vs. annular sections
    BackGround=nanmin(Chromosome.PolarContourContent);
    SumContent=nansum(Chromosome.PolarContourContent-BackGround);   
    NormContent=100*((Chromosome.PolarContourContent(Chromosome.StartLabelStopLabelBranch_CW)-BackGround)/SumContent);        
    
    %% identify which branch contains gap: this is the branch with the minimum value
    %Since we know the relative genomic locations of StartLabel, gap and StopLabel per strain, this
    %sets the orientation. Note that this depends on the experiment
    Chromosome.Orientation='Unknown';
       
    %Flip the branches if necessary  
        NormFullDistAxis=NormDistAxis;  %'historical rename'
        NormFullContent=NormContent;    
        if initval.AddFlippMarkers
            NormContent(29:31)=1.5*max(NormContent);  %marker
        end
    if (strcmp(Chromosome.Orientation,'Tails')&&initval.FlipOriention);  %incorrect orientation; flip it!
        NormFullDistAxis=fliplr(2-NormFullDistAxis);
        NormFullContent=fliplr(NormFullContent);
    end
    
    %% Sort the indices (necessary for interpolation ec.). Result is a
    %continuous axis 0....1...2 where 1 is at the StopLabel position, 0 and 2
    %depict StartLabel.
    sel=find(~isnan(NormFullDistAxis));
    NormFullContent=NormFullContent(sel); 
    NormFullDistAxis=NormFullDistAxis(sel);
     
     
    [NormFullDistAxis,idx]=sort(NormFullDistAxis);
    NormFullContent=NormFullContent(idx); 
    NormFullDistAxis=100*NormFullDistAxis;
    LC=max(NormFullDistAxis);
    [NormFullDistAxis,idx2]=unique(NormFullDistAxis);
    NormFullContent=NormFullContent(idx2);  
    
    
    
    %% 5 Up to here, all analysis is still versus the original annular axis.
    % Now, Interpolate the contents vs. outer contour length axis
    Chromosome.Ipol_OuterContourAxis=linspace(0,LC,100);
    Chromosome.Ipol_OuterContourContent=interp1(NormFullDistAxis,NormFullContent,Chromosome.Ipol_OuterContourAxis);
    
       
    %6) Section-wise content:
    %Build a BasePairAxis by integrating over the content 


     %clean a bit and re-find StopLabel position (index)
         sel=find(~isnan(NormFullContent));
         NormFullContentBP=NormFullContent(sel);    
         BasePairAxis=cumsum(NormFullContentBP);  %integrated content
         
         [BasePairAxisUnique,idx]=unique(BasePairAxis);
         BasePairAxisUnique=100*BasePairAxisUnique/max(BasePairAxisUnique);
         NormFullContentBP_Unique=NormFullContentBP(idx);
           
         %interpolate%
         PercAx=1:100;
         IpolContent=interp1(BasePairAxisUnique,NormFullContentBP_Unique,PercAx);
     
     %%       
         Chromosome.Ipol_StartLabelStopLabel_PercAxis=PercAx;
         Chromosome.Ipol_StartLabelStopLabel_Content=IpolContent;
         Chromosome.Ipol_StopLabelStartLabel_PercAxis=[];
         Chromosome.Ipol_StopLabelStartLabel_Content=[];
         