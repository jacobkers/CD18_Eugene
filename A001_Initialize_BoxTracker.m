function initval=A001_Initialize_BoxTracker(expname);

%Initialize section--------------------------------------------------------
%settings used for tracking
initval.tracklookahead=5;
initval.smoothlookahead=5;

%file handling: setup general paths (and make them if needed)
initval.data_inpath='D:\jkerssemakers\CD_Data_in\2018_Eugene\';
data_outpath='D:\jkerssemakers\Dropbox\CD_Data_out\2018_Eugene';


%define experiment-specific paths for loading and saving
switch expname
    case 'Figure2Pannel'
    initval.expi_inpath=[initval.data_inpath,'\Figure2Pannel\'];
    initval.expi_outpath=[data_outpath,'\Figure2Pannel\matlabresults\'];
    
    initval.kymofile='Kymograph_DNA.txt';
end
if ~isdir(initval.expi_outpath), mkdir(initval.expi_outpath); end

initval.roistartstop=Get_Roi_Props(expname);   %get the pre-clicked properties

         
function roistartstop=Get_Roi_Props(expname);
switch expname
    case '2019_01_22 non_interactive type'
        roistartstop(1).roino=9;
        roistartstop(1).startx=[68.155 68.343];  %x-positions
        roistartstop(1).startt=[2940 300];       %start times
        roistartstop(1).stopt=[1E6 1E6];         %stop times
        roistartstop(1).pre_t=[2840 100];        %pre-times
        roistartstop(1).trackhalfwidth=2;
        roistartstop(1).loopanalysishalfwidth=2;

        roistartstop(2).roino=31;
        roistartstop(2).startx=[65.25 67.726];
        roistartstop(2).pre_t=[332 1350];        %pre-times
        roistartstop(2).startt=[432 1450];
        roistartstop(2).stopt=[1E6 1E6];
        roistartstop(1).trackhalfwidth=2;
        roistartstop(2).loopanalysishalfwidth=2;
    case 'Figure2Pannel' 
        if 1
         roistartstop(1).roino=3;
         roistartstop(1).startx=[15.1611 72.8112];
         roistartstop(1).pre_t=[1 1];        %pre-times
         roistartstop(1).startt=[32.0724 1.9416];
         roistartstop(1).stopt=[807 158.2453];
         roistartstop(1).trackhalfwidth=[6 4];
         roistartstop(1).loopanalysishalfwidth=[6 4]; 

        else
         roistartstop(1).roino=3;
         roistartstop(1).startx=15.3287;
         roistartstop(1).pre_t=1;        %pre-times
         roistartstop(1).startt=23.7081;
         roistartstop(1).stopt=807;
         roistartstop(1).trackhalfwidth=8;
         roistartstop(1).loopanalysishalfwidth=8;

         roistartstop(2).roino=26;
         roistartstop(2).startx=51.9724;
         roistartstop(2).pre_t=200;        %pre-times
         roistartstop(2).startt=352.8061;
         roistartstop(2).stopt=5683;
         roistartstop(2).trackhalfwidth=10;
         roistartstop(2).loopanalysishalfwidth=12;

         roistartstop(3).roino=38;
         roistartstop(3).startx=57.6326;
         roistartstop(3).pre_t=1;        %pre-times
         roistartstop(3).startt=1;
         roistartstop(3).stopt=2301;
         roistartstop(3).trackhalfwidth=7;
         roistartstop(3).loopanalysishalfwidth=10;     
        end
end