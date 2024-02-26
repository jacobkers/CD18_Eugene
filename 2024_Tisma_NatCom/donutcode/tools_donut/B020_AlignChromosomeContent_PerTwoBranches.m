 function Chromosome=B020_AlignChromosomeContent_PerTwoBranches(Chromosome,StopLabel,StartLabel,initval);
% 'Use this section for a Quicksheet'
%-------------------------------------------------------------------------
    %This function returns an adapted content curve, re-interpolated to
    %align on StopLabel and StartLabel positions (if these make sense). The order of
    %contour points is set by the radial axis. We start at StartLabel
    
    % NEW fields added in this function 
    %Branches, aligned at expected label positions
    %     Ipol_StopLabelStartLabel_Content
    %     Ipol_StopLabelStartLabel_PercAxis
    %     Ipol_StartLabelStopLabel_Content
    %     Ipol_StartLabelStopLabel_PercAxis
    %     Ipol_OuterContourAxis
    %     Ipol_OuterContourContent
    %     Orientation: 'Heads'   
            %this is detrmined by the location of the
    %     glabal minimum
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
    %StartLabel-StopLabel and StopLabel-StartLabel branch 
    AngleNos=length(Chromosome.AnnularAxis);
    AngularIndices=1:AngleNos;  %angles run CCW!
    if StartLabelIndex<StopLabelIndex %zero-StartLabel-StopLabel
        Chromosome.StartLabelStopLabelBranch_CW=(StartLabelIndex:StopLabelIndex);
        Chromosome.StopLabelStartLabelBranch_CW=([StopLabelIndex:1:AngleNos 1:1:StartLabelIndex]);
    else                 %zero-StopLabel-StartLabel        
        Chromosome.StartLabelStopLabelBranch_CW=([StartLabelIndex:AngleNos 1:StopLabelIndex]);
        Chromosome.StopLabelStartLabelBranch_CW=(StopLabelIndex:1:StartLabelIndex);        
    end
   
    %% ----------------------------------------------------------------------
    %3)Get  length axes along the branches, starting out from StartLabel
    %Clockwise A
    p1x=Chromosome.CartesianContourEdge_X(Chromosome.StartLabelStopLabelBranch_CW);
    p1y=Chromosome.CartesianContourEdge_Y(Chromosome.StartLabelStopLabelBranch_CW);
    DistAxis_A=[0 [cumsum(((p1x(2:end)-p1x(1:end-1)).^2+(p1y(2:end)-p1y(1:end-1)).^2).^0.5)]'];
    NormDistAxis_A=DistAxis_A/max(DistAxis_A);    
    %Clockwise B
    p1x=Chromosome.CartesianContourEdge_X(Chromosome.StopLabelStartLabelBranch_CW);
    p1y=Chromosome.CartesianContourEdge_Y(Chromosome.StopLabelStartLabelBranch_CW);
    DistAxis_B=[0 [cumsum(((p1x(2:end)-p1x(1:end-1)).^2+(p1y(2:end)-p1y(1:end-1)).^2).^0.5)]'];
    NormDistAxis_B=DistAxis_B/max(DistAxis_B);
 
    %% ----------------------------------------------------    
    %4 normalize content of branches on whole Chromosome. Note that this is still
    %integrated vs. annular sections
    BackGround=nanmin(Chromosome.PolarContourContent);
    SumContent=nansum(Chromosome.PolarContourContent-BackGround);
    
    NormContent_A=100*((Chromosome.PolarContourContent(Chromosome.StartLabelStopLabelBranch_CW)-BackGround)/SumContent);
    NormContent_B=100*((Chromosome.PolarContourContent(Chromosome.StopLabelStartLabelBranch_CW)-BackGround)/SumContent);        
    
%% identify which branch contains gap: this is the branch with the minimum value
    %Since we know the relative genomic locations of StartLabel, gap and StopLabel per strain, this
    %sets the orientation. Note that this depends on the experiment
     Chromosome.Orientation='Unknown';
     HeadsOrTails=C010_Get_OrientationByTwoBranchGlobalMinimumLocation(NormContent_A,NormContent_B,initval);
     Chromosome.Orientation=HeadsOrTails;
       
    %Flip the branches if necessary  
        NormFullDistAxis=[NormDistAxis_A(1:end-1) 1+NormDistAxis_B];
        NormFullContent=[NormContent_A(1:end-1) NormContent_B];    
        if initval.AddFlippMarkers
            NormFullContent(29:31)=1.5*max(NormFullContent);  %marker
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
    LC=max(NormFullDistAxis);
    [NormFullDistAxis,idx2]=unique(NormFullDistAxis);
    NormFullContent=NormFullContent(idx2);  
    
    
    
    %% 5 Up to here, all analysis is still versus the original annular axis.
    % Now, Interpolate the contents vs. outer contour length axis
    Chromosome.Ipol_OuterContourAxis=linspace(0,LC,200);
    Chromosome.Ipol_OuterContourContent=interp1(NormFullDistAxis,NormFullContent,Chromosome.Ipol_OuterContourAxis);
    
       
    %6) Section-wise content:
    %Build a BasePairAxis by integrating over the content 
    %Do this for the 'A'; then 'B' branch
    %re-measure index of StopLabel

     %clean a bit and re-find StopLabel position (index)
         sel=find(~isnan(NormFullContent));
         NormFullContentBP=NormFullContent(sel);    
         BasePairAxis=cumsum(NormFullContentBP);  %integrated content
         [BasePairAxisUnique,idx]=unique(BasePairAxis);
         NormFullContentUnique=NormFullContentBP(idx);
         NormFullDistAxisUnique=NormFullDistAxis(idx);
         [~,StopLabelidx]=min(abs(NormFullDistAxisUnique-1));
         StopLabelBasePairVal=BasePairAxisUnique(StopLabelidx);
     
     
     %% Now, interpolate the A branch (passing the gap) such that it
     %adds up to the expected genomic length% . %Repeat this for the following StopLabel-StartLabel branch 
     %This amounts to force-aligning the content on the relative label
     %positions
     switch initval.straintype
         case 'type 2'
         section1=round((100-initval.StartLabelpos)+initval.StopLabelpos);  %CW1 branch length
         section2=100-section1;                         %CW2 branch length
         case 'type 1'
         section1=round(initval.StartLabelpos-initval.StopLabelpos);        %CW1 branch length
         section2=100-section1;                         %CW2 branch length
     end
         %a) StartLabel-StopLabel section
         StartLabelStopLabel_Idx=1:StopLabelidx-1;
         StartLabelStopLabel_Content=NormFullContentUnique(StartLabelStopLabel_Idx);
         OGT1=cumsum(StartLabelStopLabel_Content);
         OGT2=(OGT1-min(OGT1))/(max(OGT1)-min(OGT1))*section1;
         [OGT2,idx]=unique(OGT2);
         StartLabelStopLabel_Content=StartLabelStopLabel_Content(idx);        
         percax1=1:section1;
         OGT3=interp1(OGT2,StartLabelStopLabel_Content,percax1);
         
         %b StopLabel-StartLabel section
         StopLabelStartLabel_Idx=StopLabelidx:length(NormFullDistAxisUnique);
         StopLabelStartLabel_Content=NormFullContentUnique(StopLabelStartLabel_Idx);
         TO1=cumsum(StopLabelStartLabel_Content);
         TO2=(TO1-min(TO1))/(max(TO1)-min(TO1))*section2+section1;
         [TO2,idx]=unique(TO2);
         StopLabelStartLabel_Content=StopLabelStartLabel_Content(idx);        
         percax2=section1+1:100;
         TO3=interp1(TO2,StopLabelStartLabel_Content,percax2);
         
        
         Chromosome.Ipol_StartLabelStopLabel_PercAxis=percax1;
         Chromosome.Ipol_StartLabelStopLabel_Content=OGT3;
         Chromosome.Ipol_StopLabelStartLabel_PercAxis=percax2;
         Chromosome.Ipol_StopLabelStartLabel_Content=TO3;
       
         