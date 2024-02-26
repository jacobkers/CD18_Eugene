# cell code
----------------------------------------------
## general
* code for extended cell analysis as used in the Cees Dekker Lab, 2017-onwards
* git development [non-public] https://github.com/jacobkers/CD20_Cells
* git public: this repository
* this README largely follows the guidelines as required by Nature Communications

## public sharing
* download and unpack a .zip code of this private repository	
* remove all topics that are not relevant for the publication
* adapt this Readme file  and code accordingly
* build a .zip for Zenodo for public use

## associated publications
1. 'Direct observation of a crescent-shape chromosome in Bacillus subtilis'
Authors Miloš Tišma, Florian Patrick Bock, Jacob Kerssemakers, Aleksandre Japaridze, Stephan Gruber, Cees Dekker, accepted (2024).
2. MukBEF-dependent chromosomal organization in widened Escherichia coli
Authors Aleksandre Japaridze, Raman van Wee, Christos Gogou, Jacob WJ Kerssemakers, Cees Dekker; (2022)
3. Direct observation of independently moving replisomes in Escherichia coli
Authors Aleksandre Japaridze, Christos Gogou, Jacob WJ Kerssemakers, Huyen My Nguyen, Cees Dekker, Nat.Com (2020)
4. Direct imaging of the circular chromosome in a live bacterium
Authors Fabai Wu, Aleksandre Japaridze, Xuan Zheng, Jakub Wiktor, Jacob WJ Kerssemakers, Cees Dekker, Nat.Com (2019)

## package description
This package collects code to analyze multi-channel microscope images of cells. Patterns of interest may be extended, such as a fully labeled chromosome, or compact, as one or multiple labeled spots. The code is used to select and screen individual cells and to quantify locations of above labels. 

## package contents
* custom-written source code [Matlab] 
* test data to demonstrate the software code

## system requirements and used software
* Windows 10, 64bits 
* Matlab Version 2021 [custom source code]
* Oufti Version .... 
* Fiji(ImageJ) Version 1.52a
* Dip_image analysis software Version 2.9

## installation guide
* please follow the installation instructions provided by above package distributors
* installation of all platforms should take less than a day on a standard PC

## demo and instructions for use
* the code is divided in three main subsections (that are to be used in this order) and one extra:
	* 'crop': works on microsocopic data, outputs data organized per cells
	* 'donut': in-depth analysis on chromosome patterns and other labels ('spots')
	* 'repli': follow-up analysis on spots
	* 'topical subjects': shorter code platforms for various analyses
* see 'README.md' in each subdirectory for a detailed step-by-step descriptions, in order of appearance or importance in the analysis pipeline
* see 'README_code_annotations.md' in each subdirectory for a full collection of all custom function descriptions
* demo data for two cells is provided to illustrate workings and settings. Running this demo data typically only takes a few minutes
* for new data, follow the setup as described via the demo data 
 
## contributions
---------------------------------------------- 
* F. Wu and X. Zheng originally wrote and assembled the original code package
* J.Kerssemakers organized, re-edited and expanded the 'crop' code package
* J.Kerssemakers wrote the 'donut' and 'repli' code packages
* A. Japaridze, Raman van Wee, Christos Gogou and M.Tisma contributed to analysis design
* J.Kerssemakers and M.Tisma contributed to 'topical subjects'

## disclaimer
---------------------------------------------
This code was custom written and shared to illustrate used algorithms, analysis pathways etc. in relation to published results. The code was regularly used and tested. However, small bugs, not relevant for the data analysis as present in publications may still be around. Interested users are expected to have a sufficient knowledge of Matlab code to understand and adapt the code to their own wishes.


