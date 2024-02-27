function [initval,epth]=as_used_by_Raman_User_config_Sandro_nw_configfile(batchrunindex)
    %This function collects all local platform info, such as paths,
    %experiments and the like

    
initval.showplots=0;
initval.shortset=1000000;
initval.onlyfirstmovieframe=1;    
    
    
%% Main Paths
if ismac, initval.DirSep='/';else initval.DirSep='\';end;
initval.codepth=swap_path('F:\MukB-MatP paper\datasets\BN_CD_CellCode-master(1)\BN_CD_CellCode-master\replicode_with_filtering\');
initval.mainexperimentpath=swap_path('F:\MukB-MatP paper\datasets\Datasets_Raman\'); 
initval.processeddatapath=swap_path('F:\MukB-MatP paper\datasets\Datasets_Raman\');

initval.excelname='Overview_testruns_replidata.xlsx';
addpath(genpath(initval.codepth));


%% Experiment labels
switch batchrunindex
    %2179: ori ter
    %2816 or 18: repli (doubling rfp/cfp channels)

    case 0, initval.expi='ChristosTest';            epth=strcat('\A1-1\');
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
    case 23,    initval.expi='2019.01.17_BN2820+A22+42C+DAPI';   %5-channel Mukbef with Raman
                epth='\2820 (Wildtype, MukB-YFP)\2019.01.17 BN2820+A22+42C+DAPI\BN2820+A22+42C+DAPI.001\';
    case 24,    initval.expi='2019.02.04_BN2820+A22+DAPI001';   %5-channel Mukbef with Raman
                epth='\2820 (Wildtype, MukB-YFP)\2019.02.04 BN2820+A22+DAPI001 (good images)\';                
    case 25,    initval.expi='2019.02.04_BN2820+A22+DAPI002';   %5-channel Mukbef with Raman
                epth='\2820 (Wildtype, MukB-YFP)\2019.02.04 BN2820+A22+DAPI002\'; 
    case 26,    initval.expi='2019.02.27_BN2820+A22+DAPI';   %5-channel Mukbef with Raman
                epth='\2820 (Wildtype, MukB-YFP)\2019.02.27 BN2820+A22+DAPI\2019.02.27 BN2820+A22+DAPI\';
    case 27,    initval.expi='2019.01.17_BN2822+A22+42C+DAPI.001';   %5-channel Mukbef with Raman
                epth='\2822 (Delta MatP, MukB-YFP)\2019.01.17\BN2822+A22+42C+DAPI.001\';                
    case 28,    initval.expi='2019.02.04_BN2822+A22+DAPI002';   %5-channel Mukbef with Raman
                epth='\2822 (Delta MatP, MukB-YFP)\2019.02.04 BN2822+A22+DAPI002\';                
    case 29,    initval.expi='2019.02.27_BN2822+A22+DAPI';   %5-channel Mukbef with Raman
                epth='\2822 (Delta MatP, MukB-YFP)\2019.02.27 BN2822+A22+DAPI\2019.02.27 BN2822+A22+DAPI\';                
    case 30,    initval.expi='2019.02.27_BN2820+A22+DAPI_Tresholded_for_DNA';   %5-channel Mukbef with Raman
                epth='\2820 (Wildtype, MukB-YFP)\2019.02.27 BN2820+A22+DAPI\2019.02.27 BN2820+A22+DAPI\';                                
    case 31,    initval.expi='2019.02.27_BN2820+A22+DAPI_Tresholded_for_foci';   %5-channel Mukbef with Raman
                epth='\2820 (Wildtype, MukB-YFP)\2019.02.27 BN2820+A22+DAPI\2019.02.27 BN2820+A22+DAPI\';                                
    case 32,    initval.expi='2019.02.27_BN2822+A22+DAPI_Tresholded_for_DNA';   %5-channel Mukbef with Raman
                epth='\2822 (Delta MatP, MukB-YFP)\2019.02.27 BN2822+A22+DAPI\2019.02.27 BN2822+A22+DAPI\';                
    case 33,    initval.expi='2019.02.27_BN2822+A22+DAPI_Tresholded_for_foci';   %5-channel Mukbef with Raman
                epth='\2822 (Delta MatP, MukB-YFP)\2019.02.27 BN2822+A22+DAPI\2019.02.27 BN2822+A22+DAPI\';                
    case 34,    initval.expi='2018_10_08_MukEQ';   %5-channel Mukbef with Raman
                epth='\MukEQ\2018_10_08_MukEQ\MukEQ+DAPI.001\';

end
initval=get_exclusionlist(initval,batchrunindex);


function initval=get_exclusionlist(initval,batchrunindex);
%% exclusion list
switch batchrunindex
    case 0,    %'VersionTest';          epth='Testdata\';
        initval.exclusionlist=[{'100170'},{[]},{[]},{[]},{[]},{[]}];
    case 24,   %'Raman1, '5-channel Mukbef with Raman       
  initval.exclusionlist=[{'cell_002xy4'},{'cell_004xy4'},{'cell_003xy2'},{'cell_005xy2'},{'cell_005xy3'},{'cell_006xy2'},{'cell_007xy3'},{'cell_007xy4'},{'cell_008xy2'},{'cell_008xy4'},{'cell_009xy2'},{'cell_009xy3'}];
  case 25,   %'Raman1, '5-channel Mukbef with Raman       
  initval.exclusionlist=[{'cell_009xy1'},{'cell_009xy2'},{'cell_009xy6'},{'cell_009xy7'},{'cell_010xy2'},{'cell_010xy6'},{'cell_011xy1'},{'cell_011xy2'},{'cell_011xy4'},{'cell_012xy7'},{'cell_013xy1'},{'cell_013xy2'},{'cell_014xy4'},{'cell_015xy6'},{'cell_015xy7'},{'cell_016xy1'},{'cell_016xy3'},{'cell_017xy3'},{'cell_017xy4'},{'cell_017xy6'},{'cell_018xy1'},{'cell_018xy3'},{'cell_018xy6'},{'cell_019xy6'},{'cell_019xy7'},{'cell_020xy1'},{'cell_020xy6'},{'cell_021xy1'},{'cell_021xy8'},{'cell_022xy4'},{'cell_022xy9'},{'cell_022xy3'},{'cell_025xy9'},{'cell_026xy9'},{'cell_027xy8'}];
  case 23,   %'Raman1, '5-channel Mukbef with Raman       
  initval.exclusionlist=[{'cell_007xy7'},{'cell_008xy1'},{'cell_008xy4'},{'cell_009xy1'},{'cell_009xy8'},{'cell_010xy5'},{'cell_010xy8'},{'cell_011xy3'},{'cell_011xy6'},{'cell_011xy5'},{'cell_011xy8'},{'cell_012xy8'},{'cell_012xy9'},{'cell_013xy2'},{'cell_013xy6'},{'cell_013xy9'},{'cell_014xy3'},{'cell_014xy5'},{'cell_015xy1'},{'cell_016xy1'},{'cell_016xy2'},{'cell_016xy3'},{'cell_016xy4'},{'cell_016xy5'},{'cell_016xy6'},{'cell_016xy7'},{'cell_017xy2'},{'cell_017xy7'},{'cell_017xy8'},{'cell_017xy9'},{'cell_018xy1'},{'cell_018xy2'},{'cell_018xy4'},{'cell_018xy6'},{'cell_018xy9'},{'cell_019xy1'},{'cell_019xy6'},{'cell_020xy4'},{'cell_020xy6'},{'cell_021xy1'},{'cell_021xy4'},{'cell_021xy7'},{'cell_022xy1'},{'cell_022xy2'},{'cell_022xy6'},{'cell_023xy1'},{'cell_023xy8'},{'cell_023xy9'},{'cell_024xy2'},{'cell_024xy1'},{'cell_025xy6'},{'cell_026xy1'},{'cell_028xy7'},{'cell_028xy8'},{'cell_029xy8'},{'cell_030xy1'},{'cell_030xy5'},{'cell_030xy8'},{'cell_031xy5'},{'cell_031xy6'},{'cell_031xy7'},{'cell_032xy5'},{'cell_032xy6'},{'cell_032xy8'},{'cell_033xy7'},{'cell_034xy5'},{'cell_035xy8'},{'cell_036xy1'},{'cell_037xy5'},{'cell_037xy6'},{'cell_037xy8'},{'cell_038xy5'},{'cell_038xy6'}];
case 27,   %'Raman1, '5-channel Mukbef with Raman       
  initval.exclusionlist=[{'cell_009xy5'},{'cell_011xy1'},{'cell_011xy4'},{'cell_014xy1'},{'cell_014xy4'},{'cell_017xy3'},{'cell_019xy3'},{'cell_021xy4'}];
case 28,   %'Raman1, '5-channel Mukbef with Raman       
  initval.exclusionlist=[{'cell_007xy8'},{'cell_008xy1'},{'cell_008xy2'},{'cell_008xy4'},{'cell_009xy1'},{'cell_009xy2'},{'cell_009xy5'},{'cell_010xy1'},{'cell_010xy4'},{'cell_010xy7'},{'cell_011xy6'},{'cell_012xy3'},{'cell_013xy4'},{'cell_013xy6'},{'cell_013xy7'},{'cell_014xy1'},{'cell_014xy7'},{'cell_015xy2'},{'cell_015xy3'},{'cell_015xy7'},{'cell_016xy4'},{'cell_016xy7'},{'cell_016xy8'},{'cell_017xy4'},{'cell_020xy4'},{'cell_021xy4'}];
case 34,   %'Raman1, '5-channel Mukbef with Raman       
  initval.exclusionlist=[{'cell_002xy1'},{'cell_002xy2'},{'cell_003xy3'},{'cell_004xy2'},{'cell_005xy6'},{'cell_007xy3'},{'cell_008xy1'},{'cell_009xy1'},{'cell_009xy5'},{'cell_010xy4'},{'cell_013xy6'},{'cell_005xy1'},{'cell_005xy2'},{'cell_007xy2'},{'cell_008xy4'},{'cell_013xy3'},{'cell_014xy2'},{'cell_014xy4'},{'cell_014xy5'},{'cell_016xy4'},{'cell_016xy6'},{'cell_017xy5'},{'cell_018xy1'},{'cell_020xy4'},{'cell_020xy5'},{'cell_022xy3'},{'cell_023xy1'},{'cell_024xy1'},{'cell_025xy5'},{'cell_026xy1'},{'cell_026xy2'},{'cell_027xy2'},{'cell_026xy4'},{'cell_027xy5'},{'cell_028xy2'},{'cell_028xy5'},{'cell_029xy1'},{'cell_029xy2'},{'cell_030xy4'},{'cell_030xy2'},{'cell_030xy3'},{'cell_031xy2'},{'cell_032xy1'},{'cell_032xy4'},{'cell_034xy1'},{'cell_035xy2'},{'cell_035xy3'},{'cell_037xy1'},{'cell_037xy4'},{'cell_038xy4'},{'cell_039xy1'},{'cell_039xy4'},{'cell_041xy1'},{'cell_042xy4'},{'cell_044xy4'},{'cell_007xy7'},{'cell_008xy9'},{'cell_011xy9'},{'cell_015xy8'},{'cell_018xy8'},{'cell_019xy9'},{'cell_021xy9'},{'cell_024xy7'}];
end