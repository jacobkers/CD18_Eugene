

Code associated with the paper:

Code associated with the paper: "Bulk-surface coupling reconciles Min-protein pattern formation in vitro and in vivo"

Fridtjof Brauns1,, Grzegorz Pawlik2,, Jacob Halatek3,, Jacob Kerssemakers2, Erwin Frey1,†, Cees Dekker2,† 1Arnold Sommerfeld Center for Theoretical Physics and Center for NanoScience, Department of Physics, Ludwig-Maximilians-Universität München, Theresienstraße 37, D-80333 München, Germany 2Department of Bionanoscience, Kavli Institute of Nanoscience Delft, Delft University of Technology, Van der Maasweg 9, 2629 HZ Delft, the Netherlands 3Biological Computation Group, Microsoft Research, Cambridge CB1 2FB, UKF.B., G.P., and J.H. contributed equally to this work. †Corresponding authors. Email: frey@lmu.de, c.dekker@tudelft.nl

Code was developed in the Cees Dekker lab by Jacob Kerssemakers

Version : November_11_2020_17_58

https://zenodo.org/badge/latestdoi/312041417
__________________________________________________________
Order of listing:
main programs: project-specific shell programs
secondary programs are project-specific sub-programs
tools are small sub-functions
_________________________________________________________
       
 
Function list:
 
main analysis
________________
A001_CleanSavePatterns
A002_ResliceStacks
A003_Reslice_and_CondenseStacks
A400_Correlation
A401_Correlation_Plotter
A510_Top_and_Bottom_Correlation_XT
 
sub analysis
________________
 
tools
________________
Count_wavelength
GetBackgroundImage
JKD1_PRF_outlier_flag
MedSmooth
Scroll_ImageDirs
Split_Image
TrackXY_by_QI_Init
 
favorite tools
________________
 
Full description:
 
main analysis
________________
Name: A001_CleanSavePatterns
Parameters in: no input
Parameters out: no output
main analysis
Title: cleaning of images
Summary: this function applies background cleaning of patterns. 
 Background correction and artefact removal was done as follows: 
 first ,  for movies, frames were corrected for fluorescence bleaching 
 by normalizing each frame on its mean intensity value. 
 This corrected for the max. 20 intensity decay over long movies. 
 Next, two correction images were processed: 
 1) a ‘static background’-image Imstat was made by averaging out all 
 moving (wave pattern) features of the movie stack and removing any 
 residual background level. Thus, this image only contained static 
 fluorescent features such as specks, holes and scratches. 
 2) an ‘illumination correction’ image Imillum  was made by strongly 
 smoothening out and averaging all movie images and normalizing the result 
 to its maximum. Finally, each movie image Immovie was corrected as via the
 following image operation: Imcorrected=(Immovie- Imstat)/Imillum  .  
 This way, irregularities are suppressed and wave amplitudes on the edge 
 of each image would not be underestimated compared to the amplitudes 
 in the center of the image. 
Input: Paths are set at line ~~5. Movies should be stacked tif files
Output: cleaned tif stacks, extension 'cln'
Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
code designed & written Jacob Kerssemakers 2016 

sourcename =

    'A001_CleanSavePatterns'

Source code:A001_CleanSavePatterns line:2
__________________________________________________________ 
 
Name: A002_ResliceStacks
Parameters in: no input
Parameters out: no output
main analysis
Title: reslicing of cleaned images
Summary: this function tilts patterns: for quicker analysis, xy-t stacks
are tilted x-t and y-t, such that each saved frame is essentially a
kymograph
Input: directory with pre-cleaned movies
Output: per movie, an xt and a yt movie.
Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
code designed & written Jacob Kerssemakers 2016 

sourcename =

    'A002_ResliceStacks'

Source code:A002_ResliceStacks line:2
__________________________________________________________ 
 
Name: A003_Reslice_and_CondenseStacks
Parameters in: no input
Parameters out: no output
main analysis
Title: cleaning of images
Summary: this function tilts patterns so that each frame represents a X-T kymograph
Next, these are condensed in XY squares so the effect is, that one has a
collection of time-boxes. These can be used for top-bottom correlation
work (such as done with A510_Top_and_Bottom_Correlation_XT)
Input: 
Output:
Reference: Project: BN_CD16_Greg ; Jacob Kers 2016; Github 

sourcename =

    'A003_Reslice_and_CondenseStacks'

Source code:A003_Reslice_and_CondenseStacks line:2
__________________________________________________________ 
 
Name: A400_Correlation
Parameters in: no input
Parameters out: no output
main analysis
Title: spatial and temporal correlation analysis
Summary:  Spatial autocorrelation analysis was performed on 10 individual images per movie. 
 For each autocorrelation output image, a radial average was recorded 
 starting from the main central correlation peak. 
 The resulting spatial radial correlation curve was subjected to 
 maxima analysis. The first maximum after radius R=0  indicated 
 the most predominant distance between wave edges, 
 irrespective of propagation direction. 
 This distance was denoted as lambda. For temporal correlation, 
 we generated 20 x-t or y-t kymographs per movie 
 (10 in ‘x’ direction and 10 in ‘y’ direction) 
 evenly distributed over the middle 0.7 fraction of an image. 
 For each such kymograph, an autocorrelation analysis was performed. 
 The x=0 or y=0 line of these x-t autocorrelation maps then in effect 
 represents a temporal correlation curve averaged over all the original 
 image points on this line. Next, these correlation curves were 
 median averaged between different kymographs. 
 Thus, the final correlation curve in effect represents the average 
 temporal correlation signal sampled from 20x512 surface locations. 
 Analogous to the spatial correlation analysis, the first maximum after t=0 
 indicated a main oscillation period. 
 Experimental repeats of the same conditions (concentrations and height) 
 were median averaged.            
Input: directory with pre-cleaned movies and directory with tilted xt or yt 
versions of these movies, named following experimental
conditions. as follows:
 First letter of the filename means height, ABCDEF from lowest to highest.
 translates as heights: [2 6 8 15 25 57] microns
 second comes number which will be E/D ratio 1,2,3,4,5 , 
 translates as ratios: (0.5, 0.75. 1, 2, 3) 
 Last letter is channel (D -MinD, E -MinE)
 Second last number is experiment session
 Last number is repeat of experiment session with same conditions.
 
 Example: 'C_3_E_1_2' means: 
     height 8, ratio 1,MinE signal, session 1, repeat 2.
Output: 
 1) C_3_E_1_2_cln_overview_Height08Ratio2.00ChanERepeat01.jpg: plots
 summarizing the results of the spatial and temporal periodicity analysis
 2) A400_Correlation_results_permovie.mat:
 Data output description of Greg's 'correlation' data
   General -----------------------------------------
  MinED_mov = 1×312 struct array with fields:
      filname: name of movie
      channel: colorlabel  channel ('D' or 'E')
      height: sampe height
      EDratio: MinE-MinD ratio
      repeat: experimental repeat index of this EDratio-height combination
      xy: see below
      xt: see below
      yt: see below
      info: 
          frames2mins: 1
           pix2nm: 594
   
 Details ----------------------------------------------
 For example, MinED_mov(1).xy  contains data of the first movie.      
  Each field xy or xt contains, per movie, some details of the obained correlation curves.
 The xy correlation curves are obtained from radial averages of the
 correlation maps, the xt correlation curves are taken along the x=0 line
 (since then, we look at the autocorrelation of a wave with itself)
          suffix _mv means per movie frame (here, 10 were chosen)
          suffix _av means average of these values
 data:
          decay_val_labda_av:
          decaylength_abs_av: 
          decaylength_sign_av: 
          decayperiodicityfactor_av: 
          envelope_abs: empty
          example: [1×1 struct]: 
          mn1pos_av: position of first minimum (''anticorrelation'')
          mn1val_av: amplitude of this peak
          pk2pos_av: position of first maximum (indicative of periodicity)
          pk2val_av: amplitude of this peak
          r_axis: for plotting: axis values
          radialprofiles: for plotting: profiles
          strengthofmodulation_av: factor relating the modulation depth of the
          correlation peaks to the overall decay value, to express how
          'periodic' the waves look
 
 3) A400_Correlation_results_tabulated.mat: here, the above parameters are
 saved as columns, one row per movie
 4) A400_Correlation_results_tabulated.xlsx: as the .mat file
 5) A400_Correlation_results_matrices..jpg: visual representation of
 parameters, set in a matrix of heights and MinE/MinD ratio.
Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
code designed & written Jacob Kerssemakers 2016 

sourcename =

    'A400_Correlation'

Source code:A400_Correlation line:2
__________________________________________________________ 
 
Name: A401_Correlation_Plotter
Parameters in: no input
Parameters out: no output
main analysis
Title: matrix plotting of results
Summary: This function perfoms matrix plotting of  data processed by A400
Input: mat files in specified output dirspecified
Output: matrix plots
Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
code designed & written Jacob Kerssemakers 2016 

sourcename =

    'A401_Correlation_Plotter'

Source code:A401_Correlation_Plotter line:2
__________________________________________________________ 
 
Name: A510_Top_and_Bottom_Correlation_XT
Parameters in: no input
Parameters out: no output
main analysis
Title: load 2 stacks and correlate them image by image; allows for  timeshift
Summary: This code is used to compare MinED patterns facing each other on the top and the bottom of a sample 
Input: movies taken at different sample heights are loaded. To ave time,
these movies were tilted and resliced, so that xt-profiles can be loaded at once. 
Output: correlation data, in mat, excel and jpg
Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
code designed & written Jacob Kerssemakers 2016 

sourcename =

    'A510_Top_and_Bottom_Correlation_XT'

Source code:A510_Top_and_Bottom_Correlation_XT line:2
__________________________________________________________ 
 
 
sub analysis
________________
 
tools
________________
Name: Count_wavelength
Parameters in: resliceim
Parameters out:  Lw
tools
Title: count wave peaks in image
Input: kymograph
Output: wavelength
Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
code designed & written Jacob Kerssemakers 2019 

sourcename =

    'Count_wavelength'

Source code:Count_wavelength line:2
__________________________________________________________ 
 
Name: 
Parameters in: addvals, diagnosefilename
Parameters out: no output
tools
Title: compare different analysis routines in large analysis platforms
Summary: This function collects two datasets; for example from two
comparative analysis routines and saves the results in separate files.
to avoid passing arrays or parameters through all functions, it
initializes the container or just adds
Input: data points or arrays, plus a nr for indicating options
Output: mat files that are updated
Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
code designed & written Jacob Kerssemakers 2019 

sourcename =

  1×0 empty <a href="matlab:helpPopup char" style="font-weight:bold">char</a> array

Source code: line:2
__________________________________________________________ 
 
Name: GetBackgroundImage
Parameters in: im,backsquaregridsize
Parameters out:  im_back
tools
    Title: Get a background illumination image.
    Summary: Background Correction: take the medians of 
    sub_regions or 'tiles' and  smooth the result. The resulting image is
    taken as representative for an uneven illumination background
    Input: image, number of tiles per image width
    Output: background image, same size as original
    Reference: Cees Dekker Lab, Project: MinED; researcher Grzegorz Pawlik; 
    code designed & written Jacob Kerssemakers 2016 

sourcename =

    'GetBackgroundImage'

Source code:GetBackgroundImage line:2
__________________________________________________________ 
 
Name: JKD1_PRF_outlier_flag
Parameters in: data,tolerance,sigchange,how,sho
Parameters out: flag,cleandata
tools
    this function is meant to find a representative value for a standard
    deviation in a heavily skewed distribution (typically, flat prf with
     peaks). It calculates the standard deviation and average the prf;
     Based on these, outliers are determined and excluded for a new calculation
     of average and SD; this is repeated until sigma does not change anymore too much
     This is repeated until the new sigma does not change much anymore
    output: positions of outliers. Jacob Kers 2013 and before
    input (suggested):
         data: single data array
         tolerance: beyond how many sigmas is considered outlier (3)
         sigchange: stop iteration if relative decrease of sigma is less (0.7)
         how: consider positive outliers, or all ('all' /'positive'
         sho: (for demo only) show intermediate graphs (0)
    output: 
        flags: positions of outliers
        cleandata: data w/o outliers
      Reference: Cees Dekker Lab, 
      code designed & written Jacob Kerssemakers 2016 

sourcename =

    'JKD1_PRF_outlier_flag'

Source code:JKD1_PRF_outlier_flag line:2
__________________________________________________________ 
 
Name: MedSmooth
Parameters in: data,window,how
Parameters out:  data_smz
tools
Title: median or average smoothing
Input: array, options
Output: smoothed array
References: written by Jacob Kers, 2010 or so

sourcename =

    'MedSmooth'

Source code:MedSmooth line:2
__________________________________________________________ 
 
Name: Scroll_ImageDirs
Parameters in: impth,SearchTemplate
Parameters out:  pic_list
tools
Title: Scroll image directories
Summary:This function generates a list of full_pathfilenames and image names 
of sub-directories containing the right template string; to be used for
expanded (image) data directories
Input: source path, search template (*.[text]);
Output:
    1) structure list of pairs of full_pathfilenames and image names 
    containing the right template string
    2) list of directories that contain the specified template
References: written by Jacob Kers, 2010 or so

sourcename =

    'Scroll_ImageDirs'

Source code:Scroll_ImageDirs line:2
__________________________________________________________ 
 
Name: Split_Image
Parameters in: RawImage
Parameters out: ObjectImage, IllumImage,ImProps
tools
Title: median or average smoothing
Summary: Perform an iterative loop aimed at splitting the original image in
two images: 1) one an 'illumination field' image, that represents a 
median-smoothend, 2D curved illumination field. 2) a flattened 'objects' 
image that only deviates from zero where surface objects are present;
Input: image
Output: two images, properties
References: written by Jacob Kers, 2018 or so

sourcename =

    'Split_Image'

Source code:Split_Image line:2
__________________________________________________________ 
 
Name: TrackXY_by_QI_Init
Parameters in: firstim
Parameters out:   pretracksettings
tools
Title: setup a radial sampling map
Input: image
Output: settings and sampling coordinates
References: written by Jacob Kers, 2017 or so

sourcename =

    'TrackXY_by_QI_Init'

Source code:TrackXY_by_QI_Init line:2
__________________________________________________________ 
 
 
