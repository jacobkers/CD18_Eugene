README_details_Tisma_donutcode_March_14_2023_12_18
__________________________________________________________
This text collects lists all annotated comment sections in the code listed in this directory.
Order of listing:
main: project-specific shell programs
secondary: project-specific sub-programs
tools: small general-use sub-functions
_________________________________________________________
       
 
# Function list
 
________________
## Main analysis:
 
* A000_WF__ClickAndGo
* A010_WF_PerCell_AnalyzeAllColors
* A013_WF_PerCell_AnalyzeSpotsStandAlone
* A015_WF_PerCell_Okaying
* A0205_View_All_Channels
* A020_WF_PerCell_ChromosomeAlignment
* A040_WF_PlotGapsAndPeaks
* A050_WF_2DClusterAnalysis
* A080_WF_BuildDistanceDensityMaps
* A085_WF_BuildClusterDomainMaps
* A090_WF_SIMCorrelationDensityCurves_AbsoluteTime
* A095_WF_SIMCorrelationDensity_Relativetime
 
________________
## Sub analysis:
 
* A0001_WF_ConfigPerExperiment
* A000__WF_Get_JacobPathsandExperiments
* A030_WF_PlotSpaghettiCurves
* A055_WF_ReSampleContours
* B010_WF_BuildAveragedLabelPositions
* B050_EstimatePSF
* ChromosomeRadialSampling_WideField
* Processing_Fluorescence_SimplePatternAnalysis
 
________________
## Tools:
 
* Measure_SlopeDistances
* Get_ProfileFeatures
* Get_GuidedMax
 
# Full description
 
________________
## Main analysis:
 
Name: A000_WF__ClickAndGo
Parameters in: no input
Parameters out: no output
ClickAndGo shell
Description: Shell program for batch analysis of chromatin imaging
experiments. Check or uncheck the following sub-sections (a first
run of every step is necessary)
     * A010_WF_PerCell_AnalyzeAllColors
     * A015_WF_PerCell_Okaying(batchrunindex)
     * A020_WF_PerCell_ChromosomeAlignment(batchrunindex)
     * A030_WF_PlotSpaghettiCurves(batchrunindex)
     * A040_WF_PlotGapsAndPeaks(batchrunindex)
     * A050_WF_2DClusterAnalysis(batchrunindex)
     * A055_WF_ReSampleContours(batchrunindex)
Input: none
Output: various
References: Jacob Kerssemakers, Cees Dekker Lab, Delft
code reference:A000_WF__ClickAndGo line:2
__________________________________________________________ 
 
Name: A010_WF_PerCell_AnalyzeAllColors
Parameters in: (batchrunindex
Parameters out: no output
Title:Position and content analysis of fluorescence colour channels
General summary : The various color channels
(chromatin, c3 and c4 labeled positions) are analysed for positions and
fluorescence counts. The circular (angular) coordinate system to describe
the chromatin donut shape is set up.
Approach
Chromosome channel: Chromosome channel image is loaded. This can be the focal 
plane of a stack or its maximum projection. This image is masked with the
brightfield-based mask to remove chromosome intensity from adjacent cells. 
Next, we perform 'QI' style center-tracking on the ring structure to get 
the center and a radial map of the chromosome pattern similar to ref [1]. 
Each radial section yields a radial profile, which is subjected to peak 
analysis (position,content, width etc. 
2)c3/c4 channels are presumed to contain one spot each, these are
located following Llorente-Garcia / Reyes et al [2]. We use a squared
region, and for the 'putative spot' area pick a random
 coordinate within some pixels distance of the spot.
Input: .mat data from precursor analysis, including a mask encompassing 
the cell contour is used to exclude signal from adjacent cells.
Output: 
1)data is saved to a structure 'Chromosome' in .mat for further use; 
images are saved showing the result of the channel analysis.  
2) an Excel table listing basic geometries per cell
References
[1] M.T.J. van Loenhout, J. Kerssemakers , I. De Vlaminck, C. Dekker
Non-bias-limited tracking of spherical particles, enabling nanometer 
resolution at low magnification, Biophys. J. 102, Issue 10, 2362 (2012)
[2] Llorente-Garcia I. et al,  Biochim. Biophys. Acta. 2014;1837:811�824
code reference:A010_WF_PerCell_AnalyzeAllColors line:2
__________________________________________________________ 
 
Name: A013_WF_PerCell_AnalyzeSpotsStandAlone
Parameters in: no input
Parameters out: no output
Spot analysis per channel
Description: different color channels and the widefield channel are
allocated to the proper analysis functions, depending on the type of
pattern
Input: various
Output: various spot data .mat files
References: Jacob Kerssemakers, Cees Dekker Lab, Delft
code reference:A013_WF_PerCell_AnalyzeSpotsStandAlone line:2
__________________________________________________________ 
 
Name: A015_WF_PerCell_Okaying
Parameters in: (batchrunindex
Parameters out: no output
Title: 
Automatic multi-parameter cell screening
Summary: This program performs automatic screening on proper cell size, 
ori-ter positions etc. reports the results to the command line and saves 
a summary plot and excel table.   
Approach: the chromatin structure, label positions and cell wall positions 
are used to evaluate to what extend the cell has a well-analyzable donut. 
(in terms of size and center hole depth) and labels in a position that is roughly to be
expected based on their genomic location.
Input: .mat files from former analysis steps.
Output: 
1)a structure 'summary' containing various rejection parameters     
initialcount: counter to see how much cells are accepted  
cellsizeOK:  Cell size > 15, std <8 pixels
chromosomesizeOK: Chromosome size > 5, std <3.5 pixels
donutOK: center darker dan 0.3 of chromosome max, (0.6 for SIM)
teroriangleOK: ter-ring center-ori angle large than 60
terorilociiOK: ori and ter at least 0.3 times average
chromosome radius away from center
        totalok: Product of all filters
2)Tabular data excel, .mat and summary plots
3)A directory with only 'accepted' summary pictures of cells.
code reference:A015_WF_PerCell_Okaying line:2
__________________________________________________________ 
 
Name: A0205_View_All_Channels
Parameters in: no input
Parameters out: no output
Description: show data
input: 
    1) .mat databases from A010/013/060 replicode analysis.
    2) .mat database from a specified directory containing donut analysis
    results
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A0205_View_All_Channels line:2
__________________________________________________________ 
 
Name: A020_WF_PerCell_ChromosomeAlignment
Parameters in: (batchrunindex
Parameters out: no output
Chromosome alignment
* Summary: Mapping of chromosome and label positions and content on spatial
and genomic circular coordinates. 
Description: The initial polar coordinate analysis maps all locations 
and intensities against angle. Now, we define spatial and genomic
distances ('axes') along the semi-circular 'backbone ridge' of the chromosome and
define for both, equidistant axis points
* Approach: We start out with properties (content, width) vs. angular
positions. First, angular axes (indices) are sorted such that 1D curves
start at ori, running clockwise (optionally, we can start at ter). Next, 
the cells may be flipped or not depending in what angular order we find 
both labels and the main 'ter' gap. 
Then, this resorted angular-axis data is re-interpolated on an equidistant
distant axis that runs along the 'backbone ridge' of the chromatin
semicircle. 
Finally, the distance-mapped data is again re-interpolated, now on an axis
that has equidistant points per genomic content. For this re-interpolated
data, there is the option to correct by 'expected marker position': if for
example, the c4 'ter' label is found at 38 genomic content, clockwise
starting from ori while it is expected at 48, the genomic axis up to the
ter point is stretched accordingly. The genomic content after the ter
label (in this example, 62) is shrunk accordingly. If applied, the 
result of this procedure is that the ter label is by definition on the
expected location along the genomic axis of the chromosome.
Input: data in .mat files stored in former analysis steps. Option is to
add 'markers at fixed location along the density curves, to see where they
end up after flipping and/or alignment.
* Output: A structure 'Aligned' with subclasses 'Orig','Dist', 'BP'. Each contains
density-derived curves such as NormAxis, Density, NormDensity, 
NormCumDensity, NormCumDensityMarkercorrected.  This is saved to a .mat 
file for further analysis . In further file saving by subsequent 
analysis steps, the alignment style is included in the names
References: Jacob Kerssemakers, Cees Dekker Lab, Delft
code reference:A020_WF_PerCell_ChromosomeAlignment line:2
__________________________________________________________ 
 
Name: A040_WF_PlotGapsAndPeaks
Parameters in: (batchrunindex
Parameters out: no output
 
Peak and gap analysis of density curves
Summary: This function evaluates extrema of density curves. These are 
representative for the dense clusters and open gaps between them along the
semicircle of the chromatin pattern.
Approach: minima are detected per 1D chromatin density curves. They are 
classified as 'weak' when local density is between  50 and 100 of the
average density (minima above average are ignored). 'Strong' minima are
between 0 and 50. Similarly, maxima are detected and classified. In
addition, we evaluate how dilute or dense local chromatin is, by
asking how much distance (in spatial units) around a minimum we should travel 
along the chromatin semicircle for covering 2 genomic content. In the
absence of any corrugation, this would have to be 2. 
Input: data in .mat files stored in former analysis steps.
Output: Data is presented as scatter plots and histograms. Saved are 
Tabular excel files, .mat and summary plots.
code reference:A040_WF_PlotGapsAndPeaks line:2
__________________________________________________________ 
 
Name: A050_WF_2DClusterAnalysis
Parameters in: (batchrunindex
Parameters out: no output
Two-dimensional (image-based) cluster analysis
Summary: This function performs 2D cluster analysis on chromatin patterns of 
cells and stores the result per cell, per cluster.
Approach: we deconstruct the image by means of single-spot, psf-sized
Gaussians. This is done by subtracting such Gaussians from the image by
fist picking the brightest image point and subtracting one Gaussian of
equal peak intensity there, then repeating this action on the residual
image and continuing to do so until all intensity is covered by the subtracted 
Gaussians.
Next, the Gaussians are grouped in clusters. Each cluster consists of a
small group of Gaussians(typically 1-5) or 'components' that are within one psf 
distance of one of the others, such that the cluster forms an optically connected shape. 
Likewise, all components of one cluster are at least one psf away form
those of another.
These clusters are then parametrized by their center-of mass position,
content, radius of gyration and so on. Finally, their angular position with respect to  
 the chromatin geometrical center (as was determined via the 1D density curve analysis)
 is determined and stored.
Input: data in .mat files stored in former analysis steps.
Output: Data is presented as scatter plots and histograms. Saved are 
Tabular excel files, .mat and summary plots.
code reference:A050_WF_2DClusterAnalysis line:2
__________________________________________________________ 
 
Name: A080_WF_BuildDistanceDensityMaps
Parameters in: (batchrunindex
Parameters out: no output
Building of Dynamic 'Hi-C' analogons
Summary: This function build 'Hi-C' style 2D-plots using the circular mapping 
of the expanded chromatin 
Approach: We use the relation between spatial distance and genomic content 
along the ridge of the donut-shaped chromatin pattern (the '1D density curve').
as with an Hi-C map, the axes are genomic contact. as a measure for
contact, or mutual vicinity of genomic material, we use the negative of the
spatial distance as parameter (i.e., shorted distance means highest 'contact' signal)
and plot that in 2D.
Input: data in .mat files stored in former analysis steps.
Output: Data is saved as images. 
code reference:A080_WF_BuildDistanceDensityMaps line:2
__________________________________________________________ 
 
Name: A085_WF_BuildClusterDomainMaps
Parameters in: (batchrunindex
Parameters out: no output
Two-dimensional (image-based) cluster analysis
Summary: This function performs 2D cluster analysis on chromatin patterns of 
cells and stores the result per cell, per cluster.
Approach: we deconstruct the image by means of single-spot, psf-sized
Gaussians. This is done by subtracting such Gaussians from the image by
fist picking the brightest image point and subtracting one Gaussian of
equal peak intensity there, then repeating this action on the residual
image and continuing to do so until all intensity is covered by the subtracted 
Gaussians.
Next, the Gaussians are grouped in clusters. Each cluster consists of a
small group of Gaussians(typically 1-5) or 'components' that are within one psf 
distance of one of the others, such that the cluster forms an optically connected shape. 
Likewise, all components of one cluster are at least one psf away form
those of another.
These clusters are then parametrized by their center-of mass position,
content, radius of gyration and so on. Finally, their angular position with respect to  
 the chromatin geometrical center (as was determined via the 1D density curve analysis)
 is determined and stored.
Input: data in .mat files stored in former analysis steps.
Output: Data is presented as scatter plots and histograms. Saved are 
Tabular excel files, .mat and summary plots.
code reference:A085_WF_BuildClusterDomainMaps line:2
__________________________________________________________ 
 
Name: A090_WF_SIMCorrelationDensityCurves_AbsoluteTime
Parameters in: (batchrunindex
Parameters out: no output
Correlation analysis of movie data via 1D density curves, absolute time
Summary: We evaluate the changes in pattern of the chromatin over time, 
either by comparing 1D dennity curves with each other, or the original
images.
Approach: For each  density curve, its correlation with a reference curve  
is determined. This reference curve is chosen in two ways:
1) ''random'': randomized: the first curve of the whole dataset is taken 
as reference; irrespective of whether it belongs to the same cell or not. 
We expect such correlation to be randomly fluctuating around zero for 
most of the density curves
 2) ''movie'': the first curve of the each cell movie is taken as reference;
 Here, we expect a non-zero correlation for early frames in the movies,
 since the patterns change not entirely from movie frame to movie frame,
Both of these approaches can be done via density curves plotted against the 
spatial distance axis or the genomic axis
    
Lastly, correlation is done not via 1D density curves, but directly
from image to image. We expect here that large-scale shape changes of the
'donut' will be of more influence on the outcome.
Input: data in .mat files stored in former analysis steps.
Output: Data is presented as scatter plots and histograms. Saved are 
Tabular excel files, .mat and summary plots.
code reference:A090_WF_SIMCorrelationDensityCurves_AbsoluteTime line:2
__________________________________________________________ 
 
Name: A095_WF_SIMCorrelationDensity_Relativetime
Parameters in: (batchrunindex
Parameters out: no output
Correlation analysis of movie data via 1D density curves, relative time
Summary: We evaluate the changes in pattern of the chromatin, 
either by comparing 1D density curves with each other, or the original
images. Changes are evalauted as a function of relative time, i.e. the
frame time difference between the two curves or images.
Approach: For each frame time difference, corresponding density curves are correlated.
The results are grouped in two ways: 
 1) ''diff cell'': randomized: density cureves of different cells
 We expect such correlation to be randomly fluctuating around zero for 
 most of the density curves
 2) ''same cell'': 
 Here, we expect a non-zero correlation samll frame diffferences,
 since the patterns change not entirely from movie frame to movie frame,
 In addition, we can evaluate both type of correlation vs the 
 spatial distance axis or the genomic axis
    
Lastly, correlation is done not via 1D density curves, but directly
from image to image. We expect here that large-scale shape changes of the
'donut' will be of more influence on the outcome.
Input: data in .mat files stored in former analysis steps.
Output: Data is presented as scatter plots and histograms. Saved are 
Tabular excel files, .mat and summary plots.
code reference:A095_WF_SIMCorrelationDensity_Relativetime line:2
__________________________________________________________ 
 
 
________________
## Sub analysis:
 
Name: A0001_WF_ConfigPerExperiment
Parameters in: initval
Parameters out:  Experiment
 Configuration Per Experiment
 Experiment-specific settings to be called by 'init' file. General
 framework sets default settings, which overrides per specific experiment
 if wished.
    Alignment styles
         'c3_2Branch';  alignment of branches starts at c3 label      
         'c3_SingleBranch';  alignment of whole circle starts at c3 label
         'c4_2Branch';  alignment of bracnhes starts at c4 label
         'c4_SingleBranch';  alignment of whole circle starts at c4 label
code reference:A0001_WF_ConfigPerExperiment line:2
__________________________________________________________ 
 
Name: A000__WF_Get_JacobPathsandExperiments
Parameters in: batchrunindex
Parameters out:  initval
Settings per experiment and initialization
Summary: This function contains general and experiment-specific names and paths. 
In addition, experiment-specific stettings should be originally be  loaded from an 
'Experiment overview' excel file. To add an experiment, add entries to
this function and the excel file; use earlier entires as format example.
More recent experiments are run from matlab-listed settings in 'A0001_WF_ConfigPerExperiment'
code reference:A000__WF_Get_JacobPathsandExperiments line:2
__________________________________________________________ 
 
Name: A030_WF_PlotSpaghettiCurves
Parameters in: (batchrunindex
Parameters out: no output
Plotting of 1D density curves
Summary: This analysis produces heat maps of the collected aligned 
chromosome maps built by 'A020'and plots the result. further, the curves
are grouped in 2D 'kymographs' (position-index heat maps). 
Approach: The kymographs have the option to be sorted for various features. 
Furthermore, curves can be periodically padded to highlight features near
ori (set in the 'initval' file).
Input: pre-stored .mat data.
Output: summary plots and excel tables containing the curves.
code reference:A030_WF_PlotSpaghettiCurves line:2
__________________________________________________________ 
 
Name: A055_WF_ReSampleContours
Parameters in: (batchrunindex
Parameters out: no output
Plotting of 1D density curves
Summary: This analysis produces heat maps of the collected aligned 
chromosome maps built by 'A020'and plots the result. further, the curves
are grouped in 2D 'kymographs' (position-index heat maps). 
Approach: The kymographs have the option to be sorted for various features. 
Furthermore, curves can be periodically padded to highlight features near
ori (set in the 'initval' file).
Input: pre-stored .mat data.
Output: summary plots and excel tables containing the curves.
code reference:A055_WF_ReSampleContours line:2
__________________________________________________________ 
 
Name: B010_WF_BuildAveragedLabelPositions
Parameters in: (batchrunindex
Parameters out: no output
Title: Get and store averaged averaged label position
Summary: for movies, define an averaged albel position for a more steady
monitoring of cluster dynamics
Input: data in .mat files stored in former analysis steps.
Output: 'Cfp_Av, Rfp_Av structures
code reference:B010_WF_BuildAveragedLabelPositions line:2
__________________________________________________________ 
 
Name: B050_EstimatePSF
Parameters in: batchrunindex
Parameters out:  Psf_meas
Estimate Effective pointspread functionfor  cluster analysis
Summary: This function estimates the effective point spread function for 
a set of images.
Approach: 
Input: data in .mat files stored in former analysis steps.
Output: estimate of psf.
code reference:B050_EstimatePSF line:2
__________________________________________________________ 
 
Name: ChromosomeRadialSampling_WideField
Parameters in: im,QI, sho
Parameters out:  Chromosome
Creation of density curves by radial mapping of the chromatin pattern 
Summary: This function analyzes the semi-circular chromosome pattern via radial
mapping. From this map, 'density curves' are poduced of relative genomic
content vs. angle. Other properties, such as peak intensity and FWHM width
are stored as well.
Approach: First it performs COM centering routine on as-is chromosome pattern;
from this center, radial mapping is performed. Each radial section yields a radial
profile, which is subjected to peak analysis (position,content, width etc.
the position of the outer edge of these profiles yoields an 'edge' contour line
This center-of-mass of this contour line is then used to repeat the
whole sampling routine. This refinement is performed three times,
Input: image, quadrant interpolation settings, 
Output: structure 'Chromosome' containing density curves and other
further use
code reference:ChromosomeRadialSampling_WideField line:2
__________________________________________________________ 
 
Name: Processing_Fluorescence_SimplePatternAnalysis
Parameters in: FL
Parameters out: fluo,modelpic
Title: Processing_Fluorescence_PatternAnalysis
Approach: run it in autorun mode on provided test images to see.
Input: FL_ori: image
Output: structure 'fluo' with fields:
     area_bac: 238
     area_spot: 24
     backbone_I: [1x31 double]
     backbone_av: 4131.2
     backbone_sigma: 1007.1
     backbone_x: [1x31 double]
     backbone_y: [1x31 double]
     content_cytoplasm1: 3.5608e+005
     content_signal: 4.1605e+005
     content_spots1: 59972
     content_total: 1426534
     curve_medianofmax: [1x31 double]
     curve_medianofmax_yposses: [1x31 double]
     curve_medianofsum: [1x31 double]
     level_dark: 931.13
     level_edgetreshold: 1383.6
     level_fluotreshold: 1157.4
     level_medianofmax: 3726.6
     level_medianofsum: 13312
     level_peak: 7552
     level_sumcyto: 13312
     noise_dark: 113.11
     peak_xpos: 6
     peak_ypos: 17
     ratio_FS: 0.14414
     ratio_SN: 15.626
     wherebac: [238x1 double]
     wheredark: [631x1 double]
     wherefluo: [361x1 double]
     wherespot: [24x1 double]
 References: M. Charl Moolman*, Jacob W.J. Kerssemakers*, and Nynke H. Dekker
 Quantitative analysis of intracellular fluorescent foci in live bacteria
 Biophysical Journal, online publication September 4 (2015)
written by JacobvKerssemakers, 2012
code reference:Processing_Fluorescence_SimplePatternAnalysis line:2
__________________________________________________________ 
 
 
________________
## Tools:
 
Name: Measure_SlopeDistances
Parameters in: pic,Psf0,sho
Parameters out:  d_peakslopes
Summary: Get distances between peaks and valleys and steepest slopes
Approach: 
Input: data in .mat files stored in former analysis steps.
Output: estimate of psf.
code reference:B050_EstimatePSF line:102
__________________________________________________________ 
 
Name: Get_ProfileFeatures
Parameters in: prf,Radialweight,gap
Parameters out:  Profile
    -------------------------------------------------------------------------
    This function obtains properties of one radial profile (typically, of 256)
    first, profile is multiplied with hanning, initial max is determined there. 
    (to eliminate artefacts from neighbouring cells, especially near 
    low-intensity ter gaps etc. Then max-search is performed 
    on the original profile from this point on. A FWHM is obtained around this
    mask. Results are corrected for the radial oversampling we used and are in
    image pixel units. As a measure for local intensity, we take the sum of
    the profile but weighed per radial position to acount for the different
    number of contributing pixels (outer image pixels should count more). The
    outer edge is defined as the largest radial position with an intensity above
    35 of the peak value.
    ------------------------------------------------ Jacob Kerssemakers 2016 
code reference:ChromosomeRadialSampling_WideField line:200
__________________________________________________________ 
 
Name: Get_GuidedMax
Parameters in: Profile,MaskCenter,MaskWidth
Parameters out:  idx
    This function finds a maximum along a profile, but prefererentially in 
    a pre-specified region. This is done by iterative masking. 
code reference:Get_GuidedMax line:2
__________________________________________________________ 
 
 
# Summary
Matlab code associated with project:D:\jkerssemakers\Dropbox\CD_recent\BN_CD22_Tisma\analysis\CD20_Cells\donutcode\
Generated by Jacobs_Matlab_ManualMaker; date:March_14_2023_12_18
number of searched characters:1109739
number of searched code lines:28566
number of searched functions:495
number of included descriptions:23
 
__________________________________________________________
