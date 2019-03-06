function initval=A001_Initialize_Kymo(expname);

%Initialize section--------------------------------------------------------
%vector indicating which region-of-interest to use
initval.roilist=[3];

%settings used for tracking
initval.tracklookahead=5;
initval.smoothlookahead=5;

%file handling: setup general paths (and make them if needed)
initval.data_inpath='D:\jkerssemakers\CD_Data_in\2018_Eugene\';
initval.data_outpath='D:\jkerssemakers\Dropbox\CD_Data_out\2018_Eugene';


%define experiment-specific paths for loading and saving
switch expname
    case 'Figure2Pannel'
    initval.expi_inpath=[initval.data_inpath,'\Figure2Pannel\'];
    initval.expi_outpath=[initval.data_outpath,'\Figure2Pannel\'];
    initval.roilist=[3 26 38];
    initval.kymofile='Kymograph_DNA.txt';
end
if ~isdir(initval.expi_outpath), mkdir(initval.expi_outpath); end
