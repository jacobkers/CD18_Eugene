# crop code 
* works on microsocopic data 
* outputs data organized per cells
* output to be used by the code packages 'donut' and 'repli'
 
## general
* see the main 'README' file for installation and setup
* install dipimage (https://diplib.org/diplib-docs/dipimage_user_manual.html)

## input
* raw microscopy tiff files, see '\testdata_in\tiffs'
* deconvolution image files, see '\testdata_in\deconvolve'

## output
* a series of local output dirs
* naming as in X[020] etc. refers to the code producing it
* X020_cellcoordinate:
	celldata: cellindex, posnum ; 
	fdt1xy1 and up: celldat1, lab1
* X030_cellcrop:
	cellc1 77x116x3 data
* X040_cellgrid: summary pictures
* X040_measurement: alldata: cellc1 , cellchr, cellindex, chrindex, chrindexterincluded; oriiindex; terindex
* X050_densitydata: example: 1_cell_001t1xy1.mat contains: cellc1 82x56x3  data


## demo 
1. open X0000_AutoDataRunShell and select: user 'Jacob'; 
2. make sure the X020 to X050 are all set to '1'.
3a. set: expstotestrun=[0]:  code will process just two cells as a code check (order 1 min)
3b. set: expstotestrun=[-1]: code will process 141 cells as a data plots check (few minutes)
4. run it for a repository check


## instructions for use
1. for new image data, copy your own user file from 'X_User_config_Jacob'
2. inspect the existing case 'expstotestrun=[-1]', copy and adapt the settings accordingly
3. open X000_setpath4snapshots and add your name and your configuration file
4. open run the 'X0000_AutodataTestShell' and set it up for your user name
5. run it

## contributions
---------------------------------------------- 
* F. Wu and X. Zheng originally wrote and assembled the original code package
* J.Kerssemakers organized, re-edited and expanded the 'crop' code package

