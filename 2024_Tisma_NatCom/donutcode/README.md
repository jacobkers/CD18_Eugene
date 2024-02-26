# donut code
* works on microscopic data and crop code output 
 
## general
* see the main 'README' file for installation and setup

## input
* raw microscopy tiff files, see '\testdata_in\tiffs'
* deconvolution image files, see '\testdata_in\deconvolve'
* a series of local output dirs [X...] created by the 'crop' code 

## output
* extended mat files, excel data, jpg pictures and svg graphics. 
* in general, the prefix [A0...] refers to the code section producing it

## demo 
1. open A000_WF__ClickAndGo.
2. for a first run, make sure all A010-A055 are'set'  
3a. set: BatchrunExpArray=[0]:  code will process just two cells as a code check (order 1 min)
3b. set: BatchrunExpArray=[-1]: code will process a few tens of cells as a data plots check (few minutes)
4. run it for a repository check

## instructions for use
1. for new image data, open 'A0001_WF_ConfigPerExperiment'
2. inspect the existing examples (preferably case [-1]) and add a new experiment accordingly
3. open run the 'X0000_AutodataTestShell' and select the corresponding BatchrunExpArray.
4. run it

## contributions
* J.Kerssemakers wrote the 'donut' code package


