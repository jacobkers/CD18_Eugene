function initval=A000_ConfigExp(expi)
% 'Use this section for a Quicksheet'
    %------------------------------------------------------------------
    % This function sets user and path for the analysis programs A010/20/30   
   %------------------------------------------------------------[JK15]
 % 'End of Quicksheet section'
    
    %% some setting up
    curpth=pwd;
    cd ..; cd ..; addpath(genpath([pwd,'\common_tools\'])); cd(curpth);
    mpth='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
    %save_mpth=swap_path('Dropbox\CD_Data_out\2022_Tisma\topical_projects\2023_02_01 spot dynamics\');   
    save_mpth='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\topical_projects\2023_02_01 spot dynamics\';
    if ~isdir(save_mpth), mkdir(save_mpth); end
    
    %% paths, files
    initval.MainDataPath=[mpth '20230428_BSG4595_dynamics\'];   
    %default naming:
    initval.ImDataPath=strcat(initval.MainDataPath,expi, '_perframe\');  
    initval.SpotDataResultsName=[expi '_roi_xy.txt']; 
    initval.BleachCurveResultsName=[expi '_spots_'];
    initval.SaveDataPath=save_mpth;
    
    %% general settings
    initval.nmperpix=64;
    
    %% ROI detection
    initval.edge=20;  
    %spots closer to edge than this will not be analyzed

    initval.backsquaregridsize=100; 
    %used in background image determination

    initval.particletreshold=4;
    %number of sigmas a local maximum needs to stick out of the background 
    %to be considered

    initval.minproximity=15;
    %when spots are closer together than this, only the brightest one is kept    
    
    %% specific settings:
    switch expi
        case 'BSG4595_1s_001', initval.bleachratio_treshold=0.5;
        case 'BSG4595_1s_002', initval.bleachratio_treshold=0.5;
        case 'BSG4595_1s_004', initval.bleachratio_treshold=0.5;
        case 'BSG4595_10s_001_series1', initval.bleachratio_treshold=0.4;
        case 'BSG4595_10s_001_series2', initval.bleachratio_treshold=0.4;  
        case 'BSG4595_30s_001_series1', initval.bleachratio_treshold=0.2;
        case 'BSG4595_30s_001_series2', initval.bleachratio_treshold=0.2;
    end
    
    
    