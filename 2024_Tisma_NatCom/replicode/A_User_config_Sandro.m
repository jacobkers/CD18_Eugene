function [initval,epth]=A_User_config_Sandro(batchrunindex)
    %This function collects all local platform info, such as paths,
    %experiments and the like
    
initval.showplots=0;
initval.shortset=10E6;
initval.onlyfirstmovieframe=1;    
    
    
%% Main Paths
if ismac, initval.DirSep='/';else initval.DirSep='\';end;
initval.codepth=pwd;
initval.mainexperimentpath=swap_path('O:\MukB-MatP paper\datasets\Datasets_Raman\2820 (Wildtype, MukB-YFP)'); 
initval.processeddatapath=swap_path('O:\MukB-MatP paper\datasets\Datasets_Raman\2820 (Wildtype, MukB-YFP)\');
initval.excelname='Overview_testruns_replidata_SJ.xlsx';


 
%% Experiment labels
switch batchrunindex
    %2179: ori ter
    %2816 or 18: repli (doubling rfp/cfp channels)
    case 0,    
        initval.expi='VersionTest';          epth='Testdata\';
        initval.shortset=5;
    case 1, initval.expi='ChristosRepliSet1';       epth=strcat('\ChristosRepliSet1\');
    case 2, initval.expi='ChristosRepliSet2';       epth=strcat('\ChristosRepliSet2\');
    case 3, initval.expi='ChristosRepliSet3';       epth=strcat('\ChristosRepliSet3\');
    case 4, initval.expi='ChristosRepliSet4';       epth=strcat('\ChristosRepliSet4\');
    case 5, initval.expi='CG_23_12_2017_A1_1';      epth=strcat('\2017-12-23\A1-1\');
    case 6, initval.expi='CG_23_12_2017_A2_1';      epth=strcat('\2017-12-23\A2-1\');
    case 7, initval.expi='CG_23_12_2017_B1_1';      epth=strcat('\2017-12-23\B1-1\');
    case 8, initval.expi='CG_23_12_2017_B3_1';      epth=strcat('\2017-12-23\B3-1\');
    case 9, initval.expi='CG_23_12_2017_B3_2';      epth=strcat('\2017-12-23\B3-2\');
    case 10,initval.expi='2179-18-10-2017_crop1';   epth=strcat('\20180524_Bug\2179-18-10-2017_crop1\'); 
    case 11,initval.expi='2017-10-18-A1-2';         epth=strcat('\2017-10-18\A1-2\');
        %Replipaper 2019:
    case 12,    initval.expi='2018-04-09-2816-long cells_crop'; 
                epth=strcat('\2018-04-09-2816-long cells_crop\'); 
    case 13,    initval.expi='2018-05-23-2818+A22crop'; 
                epth=strcat('\2018-05-23-2818+A22\2818-A22_crop\');
    case 14,    initval.expi='2018-05-23-2818+A22full'; 
                epth=strcat('\2018-05-23-2818+A22 full\2818-A22\');
    case 15,    initval.expi='2018_05_28 2818_woA22_longcells'; 
                epth=strcat('\2018_05_28 2818_woA22_longcells\');
    case 16,    initval.expi='2018-05-29-2818-30C-short'; 
                epth=strcat('\2018-05-29-2818-30C-short\');
    case 17,    initval.expi='2018-05-28-2818-42C-short_t1'; 
                epth=strcat('\2018-05-28-2818-42C-short\t1\');
    case 18,    initval.expi='2018-05-28-2818-42C-short_t13\'; 
                epth=strcat('\2018-05-28-2818-42C-short\t13\');
    case 19,    initval.expi='2018-05-28-2818-42C-short_t25'; 
                epth=strcat('\2018-05-28-2818-42C-short\t25\');
    case 20,    initval.expi='2018-05-23-2818+A22 prescreened'; 
                epth=strcat('\2018-05-23-2818+A22 prescreened\');
    case 22,    initval.expi='ChannelSwapTest';   %5-channel Mukbef with Raman
                epth='\My_5ChanTestdata\';
    case 23,    initval.expi='Test';   %5-channel Mukbef with Raman
                epth='\2820\2018-06-14-2820-DAPI\2018-06-14-2820-DAPI\';
end
