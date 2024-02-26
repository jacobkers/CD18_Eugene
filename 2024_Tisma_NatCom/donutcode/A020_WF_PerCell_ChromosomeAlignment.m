function A020_WF_PerCell_ChromosomeAlignment(batchrunindex)
%JWJK_A:-------------------------------------------------------------------
%Chromosome alignment
%* Summary: Mapping of chromosome and label positions and content on spatial
%and genomic circular coordinates. 
%Description: The initial polar coordinate analysis maps all locations 
%and intensities against angle. Now, we define spatial and genomic
%distances ('axes') along the semi-circular 'backbone ridge' of the chromosome and
%define for both, equidistant axis points
%* Approach: We start out with properties (content, width) vs. angular
%positions. First, angular axes (indices) are sorted such that 1D curves
%start at ori, running clockwise (optionally, we can start at ter). Next, 
%the cells may be flipped or not depending in what angular order we find 
%both labels and the main 'ter' gap. 
%Then, this resorted angular-axis data is re-interpolated on an equidistant
%distant axis that runs along the 'backbone ridge' of the chromatin
%semicircle. 
%Finally, the distance-mapped data is again re-interpolated, now on an axis
%that has equidistant points per genomic content. For this re-interpolated
%data, there is the option to correct by 'expected marker position': if for
%example, the c4 'ter' label is found at 38% genomic content, clockwise
%starting from ori while it is expected at 48%, the genomic axis up to the
%ter point is stretched accordingly. The genomic content after the ter
%label (in this example, 62%) is shrunk accordingly. If applied, the 
%result of this procedure is that the ter label is by definition on the
%expected location along the genomic axis of the chromosome.
%Input: data in .mat files stored in former analysis steps. Option is to
%add 'markers at fixed location along the density curves, to see where they
%end up after flipping and/or alignment.
%* Output: A structure 'Aligned' with subclasses 'Orig','Dist', 'BP'. Each contains
%density-derived curves such as NormAxis, Density, NormDensity, 
%NormCumDensity, NormCumDensityMarkercorrected.  %This is saved to a .mat 
%file for further analysis . In further file saving by subsequent 
%analysis steps, the alignment style is included in the names
%References: Jacob Kerssemakers, Cees Dekker Lab, Delft
%:JWJK_A-------------------------------------------------------------------

if nargin<1,batchrunindex=11.1;end

%B010_WF_BuildAveragedLabelPositions(batchrunindex);

actions.showit=0;

initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex);

FewCells=0;  %Demo purposes
if FewCells
    %initval.Cell_Labels=[{'1000002_t04'}; {'1000002_t05'}];  %11.1
     initval.Cell_Labels=[{'200020'}; {'200022'}];           %62
     initval.Cell_Labels=[{'100578'}; {'100578'}];           %9
end

allframes=length(initval.Cell_Labels);
close all;
pix2nm=64;

%for jj=1:allframes
 for jj=1:allframes
    cellno=char(initval.Cell_Labels{jj});
    
    CellName=strcat('ResultsOfCell',cellno,'.mat');
    disp(strcat(CellName,'-Aligning genomic and distance axes..',num2str(allframes-jj+1),'cells to go'));
    MatFilePath=strcat(initval.resultpath,'ResultsPerCellMatlab',initval.DirSep);
    load(strcat(MatFilePath,CellName));  
    if GeneralCellProps.Okayproduct|initval.PassAllCells
        switch initval.alignmodus
            case 'c3_2Branch',     StartLabel=c3; StopLabel=c4; MarkerLabel=c4;
            case 'c4_2Branch',     StartLabel=c4; StopLabel=c3; MarkerLabel=c3;          
            case 'c3_SingleBranch',StartLabel=c3; StopLabel=c3; MarkerLabel=c4;
            case 'c4_SingleBranch',StartLabel=c4; StopLabel=c4; MarkerLabel=c3;
            case 'c3Av_SingleBranch',
                StartLabel=c3_AvProps; 
                StopLabel=c3_AvProps; 
                MarkerLabel=c4_AvProps;
           case 'c3Av_2Branch',
                StartLabel=c3_AvProps; 
                StopLabel=c4_AvProps; 
                MarkerLabel=c4_AvProps;
        end
        
        %% Perform alignment in three steps:
        %1) shuffle indices to get branches according to ori, ter position,
        %flip if option set. results stored in sub-field 'Ori'
        %2) Re-sample properties of interest to an equidistant axis;
        %calculate the distance along which the ter-marker is found;
        %Use the 'ori field alone        
        Aligned=EQUIANGLE_ORI_OrganizeBranchOrder(Chromosome,StopLabel,StartLabel,MarkerLabel,initval);                   
        Aligned=EQUIDIST_ResampleAlongRidge(Aligned,CellName);
        Aligned=EQUIBASEPAIR_ResampleAlongRidge(Aligned,initval); 
        %Aligned=orderfields(Aligned);
        save(strcat(MatFilePath,CellName),'Aligned', '-append');
        dum=1;
    end
end
disp('done');

function Chromosome=Add_FlipMarkers(Aligned,Chromosome) 
        LP=length(Aligned.Orig.AllIndices);
        MarkIdx=Aligned.Orig.AllIndices(ceil(LP/4));  %take the index at one quarter from the chosen label start position
        lo=max([MarkIdx-2 1]); hi=min([MarkIdx+2 LP]);
        Chromosome.PolarContourPeakVal(lo:hi)=max(Chromosome.PolarContourPeakVal);
 
        
function Aligned=EQUIANGLE_ORI_OrganizeBranchOrder(Chromosome,StopLabel,StartLabel,MarkerLabel,initval)
      % 'Use this section for a Quicksheet'
%-------------------------------------------------------------------------  
    %function will perform sorting indices according to branch A and B,
        %clockwise; will flip if necesarry; will add 'flip markers' if necessary 
   %------------------------------------------------
% 'End of Quicksheet section'
%JacobKers2017
Aligned=struct('Orig',[]);

%% Check number of label positions to use: one or two 
if StartLabel.spotPosAngle==StopLabel.spotPosAngle,
    Aligned.OneBranch=1;
else
    Aligned.OneBranch=0;
end
%1) find start(StartLabel) and stop(StopLabel) points on contour curve
    [minStopLabeldist,StopLabelIndex]=min...   %annular index of nearest point
    (((Chromosome.CartesianContourMax_X-StopLabel.spotX).^2+...
    (Chromosome.CartesianContourMax_Y-StopLabel.spotY).^2).^0.5);

    [minStartLabeldist,StartLabelIndex]=min...   %annular index of nearest point
    (((Chromosome.CartesianContourMax_X-StartLabel.spotX).^2+...
    (Chromosome.CartesianContourMax_Y-StartLabel.spotY).^2).^0.5);

    [minSMarkerLabeldist,MarkerLabelIndex]=min...   %annular index of nearest point
    (((Chromosome.CartesianContourMax_X-MarkerLabel.spotX).^2+...
    (Chromosome.CartesianContourMax_Y-MarkerLabel.spotY).^2).^0.5);
    
    Aligned.Orig.StartLabelIndex=StartLabelIndex;
    Aligned.Orig.StopLabelIndex=StopLabelIndex;    
    Aligned.Orig.MarkerLabel=MarkerLabel;
    Aligned.Orig.MarkerLabel.OriIndex=MarkerLabelIndex;

    %% 2) find angular indices in the right order to define an always clockwise
    %StartLabel-StopLabel and StopLabel-StartLabel branch (the latter maybe empty) 
    AngleNos=length(Chromosome.AnnularAxis);
    AngularIndices=1:AngleNos;  %angles run CCW!
    if ~Aligned.OneBranch
        if StartLabelIndex<StopLabelIndex %zero-StartLabel-StopLabel
            Aligned.Orig.StartStopIndices=(StartLabelIndex:StopLabelIndex);
            Aligned.Orig.StopStartIndices=([StopLabelIndex:1:AngleNos 1:1:StartLabelIndex]);
            %Note that overlapping start and end points points are removed:
            Aligned.Orig.AllIndices=[Aligned.Orig.StartStopIndices Aligned.Orig.StopStartIndices(2:end-1)];
        end
        if StartLabelIndex>=StopLabelIndex  %zero-StopLabel-StartLabel 
            Aligned.Orig.StartStopIndices=([StartLabelIndex:AngleNos 1:StopLabelIndex]);   
            Aligned.Orig.StopStartIndices=(StopLabelIndex:1:StartLabelIndex);
            %Note that overlapping start and end points points are removed:
            Aligned.Orig.AllIndices=[Aligned.Orig.StartStopIndices Aligned.Orig.StopStartIndices(2:end-1)];
        end
    else       %SINGLE BRANCH         
            Aligned.Orig.StartStopIndices=[[StartLabelIndex:AngleNos] [1:StopLabelIndex-1]];           
            Aligned.Orig.StopStartIndices=[];  
            Aligned.Orig.AllIndices=Aligned.Orig.StartStopIndices;
    end
    
     %% 2a Get orientation and optionally flip the indices accordingly:
    %identify which branch contains gap: this is the branch with the minimum value
    %Since we know the relative genomic locations of StartLabel, gap and StopLabel per strain, this
    %sets the orientation. Note that this depends on the experiment. Indices are shuffled optionally and accordingly 
    
    if initval.AddFlippMarkers,Chromosome= Add_FlipMarkers(Aligned,Chromosome); end   
    
    if ~Aligned.OneBranch        
        Aligned=C010_Get_OrientationByTwoBranchGlobalMinimumLocation(Chromosome,Aligned,initval);
    else  %single branch alignment according to different protocols
        switch initval.FlipModus, 
            case 'UseGlobalMin' 
                    Aligned=C015_Get_OrientationBySingleBranchGlobalMinimumLocation(Chromosome,Aligned,initval);           
        end
    end
    
    
    %1)Length along contour: Get  length along the outer contour, starting out from StartLabel;      
    p1x=Chromosome.CartesianContourMax_X(Aligned.Orig.AllIndices);
    p1y=Chromosome.CartesianContourMax_Y(Aligned.Orig.AllIndices);
    
    DistAxis=[0 [cumsum(((p1x(2:end)-p1x(1:end-1)).^2+(p1y(2:end)-p1y(1:end-1)).^2).^0.5)]'];
    NormDistAxis=100*DistAxis/nanmax(DistAxis);
    
    %Density=Chromosome.PolarContourPeakVal(Aligned.Orig.AllIndices);
    Density=Chromosome.PolarContourContent(Aligned.Orig.AllIndices);
    
    %build a single-pass property matrix
    %GeneralPropertyArray: rows: content, FWHM, ridgelineradius,peakvalue
    GeneralPropertyArray=[
        Chromosome.PolarContourContent(Aligned.Orig.AllIndices);    
        Chromosome.PolarContourFWHM(Aligned.Orig.AllIndices);...
        (Chromosome.PolarContourMax(Aligned.Orig.AllIndices))';...
        Chromosome.PolarContourPeakVal(Aligned.Orig.AllIndices)];
        
    
    
    CumDensity=cumsum(Density); CumDensity=CumDensity-CumDensity(1);
    NormCumDensity=100*CumDensity/nanmax(CumDensity);
    
    Aligned.Orig.AngleAsUsed=Chromosome.AnnularAxis(Aligned.Orig.AllIndices);
    Aligned.Orig.DistVsAngleAsUsed=DistAxis;
    Aligned.Orig.NormDist=NormDistAxis;
    Aligned.Orig.NormCumDensity=NormCumDensity;
    Aligned.Orig.Density=Density;
    Aligned.Orig.GeneralPropertyArray=GeneralPropertyArray;
    
    
    %Now, find the distence and content the marker is associated with
    curidx=find(Aligned.Orig.AllIndices==Aligned.Orig.MarkerLabel.OriIndex);
    curidx=max(curidx);  %patch for 0/100 entry
    %this is the location in the re-shuffledindex-array
    Aligned.Orig.MarkerLabel.MarkerDistPix=DistAxis(curidx);
    Aligned.Orig.MarkerLabel.MarkerDistPerc=NormDistAxis(curidx);
    Aligned.Orig.MarkerLabel.MarkerBPPerc=NormCumDensity(curidx);   
    switch initval.straintype
       case 'type 2', perc=round((100-initval.StartLabelpos)+initval.StopLabelpos);  %CW1 branch length
       case 'type 1', perc=round(initval.StartLabelpos-initval.StopLabelpos);        %CW1 branch length
    end
    Aligned.Orig.MarkerLabel.MarkerBPPercExpected=perc;
    
    
    
 

 function Aligned=EQUIDIST_ResampleAlongRidge(Aligned,CellName);
% 'Use this section for a Quicksheet'
%-------------------------------------------------------------------------
    %Map properties to distance axis
    %------------------------------------------------
% 'End of Quicksheet section'properties to 
%JacobKers2016
    
   %Obtain curves vs equidistant contour points, normalized to 100
    [Aligned,ResampledDensity]=C025_FromEquiAngle2EquiDist(Aligned,Aligned.Orig.Density,CellName);    
    Aligned.Dist.Density=ResampledDensity;
    Aligned.Dist.NormDensity=ResampledDensity/nansum(ResampledDensity)*100;
   
    [Aligned,ResampledNormCumDensity]=C025_FromEquiAngle2EquiDist(Aligned,Aligned.Orig.NormCumDensity,CellName);
    Aligned.Dist.NormCumDensity=ResampledNormCumDensity;   
    
    [Aligned,ResampledGeneralPropertyArray]=C025_FromEquiAngle2EquiDist(Aligned,Aligned.Orig.GeneralPropertyArray,CellName);    
    Aligned.Dist.GeneralPropertyArray=ResampledGeneralPropertyArray;
    
    
   function Aligned=EQUIBASEPAIR_ResampleAlongRidge(Aligned,initval);         
    %% Build Basepairaxis-content; starting from interpolated distance-based curves
    %marker-pos correction is done here.
   [Aligned,ResampledDensity]=C035_FromEquiDist2EquiContent(Aligned,Aligned.Dist.Density,initval);
    Aligned.BP.NormDensity=ResampledDensity/nansum(ResampledDensity)*100;
    Aligned.BP.Density=ResampledDensity;
     
    [Aligned,ResampledGeneralPropertyArray]=C035_FromEquiDist2EquiContent(Aligned,Aligned.Dist.GeneralPropertyArray,initval);
    Aligned.BP.GeneralPropertyArray=ResampledGeneralPropertyArray;
    
 
    
     
 
   
    
         





