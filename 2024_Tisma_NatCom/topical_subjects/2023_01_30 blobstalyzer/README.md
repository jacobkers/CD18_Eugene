# blobstalyzer
----------------------------------------------
## general
* code for user-dependent blind classification of cell patterns
* git development [non-public] https://github.com/jacobkers/CD20_Cells
* this README largely follows the guidelines as required by Nature Communications


## associated publication
1. 'Direct observation of a crescent-shape chromosome in Bacillus subtilis'
Authors Miloš Tišma, Florian Patrick Bock, Jacob Kerssemakers, Aleksandre Japaridze, Stephan Gruber, Cees Dekker, submitted (2023).


## package description
The program detects separate chromosome patterns from four different data sets and identifies single-cell regions-of-interest [ROIs].
Then, it  offers just one, randomly picked from any of the input images to a user. The user does not know which experiment the ROI stems from and is asked to classify it among a number of shape classes ('donut' etc. This is repeated a fixed number of times, typically 100. The result is a distribution histogram per experiment, per type of shape. The cycle is repeated 10x to increase statistics and evaluate single-user variability. Next, this whole sequence can be duplicated by a second user, or more. Users do not communicate on judgment before performing the task. Thus, any bias depending on the judgment of a single user is shown by the statistically significant differences in distributions

## package contents
* custom-written source code [Matlab] 
* test data to demonstrate the software code

## system requirements and used software
* Windows 10, 64bits 
* Matlab Version 2021 [custom source code]

## installation guide
* please follow the installation instructions provided by above package distributors
* installation of all platforms should take less than a day on a standard PC

## demo and instructions for use
1. run 'Blobstalyzer' and follow the instructions
2. repeat for as many users as wished
3. run 'collect_blobstalysis' (adapted to proper user names) to collect the results
 
## contributions
---------------------------------------------- 
* J.Kerssemakers wrote code and performed analysis
* M. Tisma collected data and performed analysis

## disclaimer
---------------------------------------------
This code was custom written and shared to illustrate used algorithms, analysis pathways etc. in relation to published results. The code was regularly used and tested. However, small bugs, not relevant for the data analysis as present in publications may still be around. Interested users are expected to have a sufficient knowledge of Matlab code to understand and adapt the code to their own wishes.




