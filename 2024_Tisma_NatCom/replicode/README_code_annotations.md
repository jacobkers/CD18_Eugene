README_details_Tisma_replicode_March_14_2023_12_18
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
 
* A000_LazyShell
* A000_Repli_Init
* A001_Build_MovieList
* A001_Erase_and_Build_ResultDirs
* A005_Build_MovieList
* A010_WF_PerCell_AnalyzeCellShapeStandAlone
* A013_WF_PerCell_AnalyzeSpots
* A055_WF_2DClusterAnalysis_Standalone
* A060_Colocalization
* A062_MukB_Dna_overlap
* A070_Shell_Collect_results
 
________________
## Sub analysis:
 
* Form_LookingBackPairs
* Form_LookingBackPairs
* Get_labelpair_distances
* Get_ori_ori_ter_distances
* Processing_Fluorescence_SimplePatternAnalysis
 
________________
## Tools:
 
* PeelblobsFromImage
 
# Full description
 
________________
## Main analysis:
 
Name: A000_LazyShell
Parameters in: no input
Parameters out: no output
Shell 
Description: Shell program to run sequential steps of main code for
replication analysis. Runs all main programs in set order; former runs can be
re-used by setting the switches 'if 1'
input: prior 'crop code' analysis data.
output:  new .mat database is generated and expanded per step.
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A000_LazyShell line:2
__________________________________________________________ 
 
Name: A000_Repli_Init
Parameters in: batchrunindex,usr
Parameters out:  initval
Initialization
Description: set main paths of data to be analyzed. Settings per
experiment and initialization 
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A000_Repli_Init line:2
__________________________________________________________ 
 
Name: A001_Build_MovieList
Parameters in: no input
Parameters out: no output
Build movielist
Description: identifies sets of frames per movie based on the filenames
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A001_Build_MovieList line:2
__________________________________________________________ 
 
Name: A001_Erase_and_Build_ResultDirs
Parameters in: (initval
Parameters out: no output
Description: For first time use of experiment; Check well before you use it!!
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A001_Erase_and_Build_ResultDirs line:2
__________________________________________________________ 
 
Name: A005_Build_MovieList
Parameters in: (batchrunindex
Parameters out: no output
Build movielist
Description: identifies sets of frames per movie based on the filenames
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A005_Build_MovieList line:2
__________________________________________________________ 
 
Name: A010_WF_PerCell_AnalyzeCellShapeStandAlone
Parameters in: (initval
Parameters out: no output
Description: obtains general cell shape data
input: .mat database from crop code analysis.
output:.mat database is generated and expanded at later steps.
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A010_WF_PerCell_AnalyzeCellShapeStandAlone line:2
__________________________________________________________ 
 
Name: A013_WF_PerCell_AnalyzeSpots
Parameters in: (initval
Parameters out: no output
Description: analyze local lables (typically ori and ter) in their
respective channels
input: .mat database from A010 analysis.
output:.mat database is expanded, cell images are saved; some tables. All
with 'A013' label.
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A013_WF_PerCell_AnalyzeSpots line:2
__________________________________________________________ 
 
Name: A055_WF_2DClusterAnalysis_Standalone
Parameters in: (initval
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
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A055_WF_2DClusterAnalysis_Standalone line:2
__________________________________________________________ 
 
Name: A060_Colocalization
Parameters in: (initval
Parameters out: no output
Description: compare channels pair-wise and determine degree of overlap of signal 
input: .mat database from A010/013 analysis.
output: .mat database labeled A060
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A060_Colocalization line:2
__________________________________________________________ 
 
Name: A062_MukB_Dna_overlap
Parameters in: (initval
Parameters out: no output
Description: process colocalization data 
input: 
    1) .mat databases from A010/013/060 replicode analysis.
    2) .mat database from a specified directory containing donut analysis
    results
output: directory labeled A062 with images, table, database
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A062_MukB_Dna_overlap line:2
__________________________________________________________ 
 
Name: A070_Shell_Collect_results
Parameters in: no input
Parameters out: no output
Shell 
Description: Shell program to run sequential steps of main code for
replication analysis. Runs all main programs in set order; former runs can be
re-used by setting the switches 'if 1'
input: prior 'crop code' analysis data.
output:  new .mat database is generated and expanded per step.
Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
code reference:A070_Shell_Collect_results line:2
__________________________________________________________ 
 
 
________________
## Sub analysis:
 
Name: Form_LookingBackPairs
Parameters in: Rfp,Rfp_former,Cell,Cell_former,options
Parameters out: ori_pairsX,ori_pairsY,FormerNearXYCoords
    Title: Form time-pairs
    Summary: this function forms coordinate pairs of a spot and its nearest spot in the
    former image; it is assumed these represent the same object in time     
    ''Rfp'' with fields: spotY, spotX,spotContent, same for fomrer image    
    Output: separate X and Y coordinates of associated pairs and an xy
    list for future use
code reference:P100_AnalyzeSpotGeometries line:64
__________________________________________________________ 
 
Name: Form_LookingBackPairs
Parameters in: Rfp,Rfp_former,Cell,Cell_former,options
Parameters out: ori_pairsX,ori_pairsY,FormerNearXYCoords
    Title: Form time-pairs
    Summary: this function forms coordinate pairs of a spot and its nearest spot in the
    former image; it is assumed these represent the same object in time     
    ''Rfp'' with fields: spotY, spotX,spotContent, same for fomrer image    
    Output: separate X and Y coordinates of associated pairs and an xy
    list for future use
code reference:A100_AnalyzeSpotGeometries line:64
__________________________________________________________ 
 
Name: Get_labelpair_distances
Parameters in: Rfp,Cell,options
Parameters out:  pair_table
    Title:   Classify replication
    Summary: A single cell image is judged by various distance parameters  
    Approach: each ori spot is associated with its nearest-neighbour in
    the same frame. Their common COM point is determined. this point is
    related to the center-of mass of the cell, and/or the ter position
    Input: 
    Cfp:   Rfp with fields: spotY, spotX,spotContent,
    Cell:  Area:Centroid MajorAxisLength MinorAxisLength Orientation BW Edge 
    Output: 
        'dist' table with typical pair distances, 
        'pair' table with coordinates of nearest-neighbour  pairs
        [ii xo xo2 xp yo yo2 yp xm ym d_oo]
        In both cases, indices run by those of the ori spots
code reference:Get_labelpair_distances line:2
__________________________________________________________ 
 
Name: Get_ori_ori_ter_distances
Parameters in: Cfp,Rfp,Cell
Parameters out:  ori_ori_ter_table
    Title:   Classify replication
    Summary: A single cell image is judged by various distance parameters  
    Approach: each ori spot is associated with its nearest-neighbour in
    the same frame. Their common COM point is determined. this point is
    related to the center-of mass of the cell, and/or the ter position
    Input: 
    Cfp:   Rfp with fields: spotY, spotX,spotContent,
    Cell:  Area:Centroid MajorAxisLength MinorAxisLength Orientation BW Edge 
    Output: 
        'dist' table with typical pair distances, 
        'pair' table with coordinates of nearest-neighbour  pairs
        [ii xo xo2 xp yo yo2 yp xm ym]
        In both cases, indices run by those of the ori spots
code reference:Get_ori_ori_ter_distances line:2
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
 
Name: 
Parameters in: no input
Parameters out: no output
    Alignment styles
         'Rfp_2Branch';  alignment of branches starts at rfp label      
         'Rfp_SingleBranch';  alignment of whole circle starts at rfp label
         'Cfp_2Branch';  alignment of bracnhes starts at cfp label
         'Cfp_SingleBranch';  alignment of whole circle starts at cfp label
code reference:A0001_WF_ConfigPerExperiment line:3
__________________________________________________________ 
 
Name: PeelblobsFromImage
Parameters in: im,Psf,ChipFract,RelChangeStop,sho,SepSigs
Parameters out:  AllSpotProps
Title: Multi-peak fitter
Summary: This function subtracts 2D Gaussians of fixed width from an image until a stop criterion
is met. the gaussians can then be used for reconstructing the curve in
main components, for example defined by optical separation.
Output: AllSpotProps=[spotcount Peak Xpos Ypos Psf ThisSpotFraction(spotcount) CoveredFraction(spotcount) RelChange]];  
Project: BN_CD16_Greg, Fabai ; JacobKers 2016
code reference:PeelblobsFromImage line:2
__________________________________________________________ 
 
 
# Summary
Matlab code associated with project:D:\jkerssemakers\Dropbox\CD_recent\BN_CD22_Tisma\analysis\CD20_Cells\replicode\
Generated by Jacobs_Matlab_ManualMaker; date:March_14_2023_12_18
number of searched characters:558987
number of searched code lines:13548
number of searched functions:255
number of included descriptions:18
 
__________________________________________________________
