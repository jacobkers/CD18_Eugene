function initval=A_User_config_Jacob(batchrunindex)
%Description: %This function collects all local platform info, such as paths,
%experiments and the like
%input: index referring to experiment
%output:general and specific intialization
%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20

%% RUNTIME settings:
initval.showplots=0;
initval.shortset=20;
initval.onlyfirstmovieframe=0;      

%% DEFAULT per exp (can be overwritten below or when loading excel)
initval.searchlabel='c';
initval.cellmaskname='cellmask';
initval.analyze_ori_ter=1; 
initval.analyze_replisome=0; 
initval.Psf_est=2.9; 
initval.UseMeasuredPSFforClusterAnalysis=0;
initval.cluster_sep_sigs=2; 
initval.nmperpixel=65; 
initval.orimaskradius=10;
initval.axislimit=5;
%channel info:
initval.chan_ID=[1 0 2 3 4];  %'NEW' treat as phase msk chr ori ter
                 
initval.tresholdperchannel=[NaN 0.3  0.7 0.7 0.4]; 
%(comment out to use auto) [WF DNA   Ori Ter MukB]

initval.channelorder=2;
initval.numberofchannels=5;  
initval.pth_excelfile='none';

%check where to find data (in order of acces speed)
drivepath=set_drivepath;

%% PER EXPERIMENT
switch batchrunindex  
    case -3
            initval.expi='Tisma_4595_IPTG_2h_SyG_test';
            initval.pth_crop=swap_path('CD_Data_in\2022_Tisma\20220421_4595_tests_analysis\tiffcz_4595_IPTG_2h_SyG\');               
            initval.pth_donut=swap_path('Dropbox\CD_Data_out\2022_Tisma\BatchanalysisResults\Results_Tisma_4595_IPTG_2h_SyG_test_Rfp_SingleBranch_NoFlip\');                       
            initval.pth_repli=swap_path('Dropbox\CD_Data_out\2022_Tisma\\BatchanalysisResults_Repli\Tisma_4595_IPTG_2h_SyG_test\');         
            initval.chan_ID=[1 0 2 3 4];  %'treat as phase msk chr ori ter
    case -2
            initval.expi= 'Tisma_4595_IPTG_2h_test'  ;   
            initval.pth_crop=swap_path('CD_Data_in\2022_Tisma\20220421_4595_tests_analysis\tiffcz_4595_IPTG_2h\');               
            initval.pth_donut=swap_path('Dropbox\CD_Data_out\2022_Tisma\BatchanalysisResults\Results_Tisma_4595_IPTG_2h_test_Rfp_SingleBranch_NoFlip\');                       
            initval.pth_repli=swap_path('Dropbox\CD_Data_out\2022_Tisma\BatchanalysisResults_Repli\Tisma_4595_IPTG_2h_test\');         
            initval.chan_ID=[1 0 2 3 4];  %'treat as phase msk chr ori ter       
    case -1,initval.expi='repo_test';
            codepth=pwd; 
            cd .. ; 
            initval.pth_crop= [pwd, '\testdata_in\'];
            initval.pth_donut=[pwd, '\testdata_out\BatchanalysisResults\Results_repo_test_Rfp_2Branch_NoFlip\']; 
            initval.pth_repli=[pwd, '\testdata_out\BatchanalysisResults_Repli\repo_test\'];
            cd(codepth); 
            initval.numberofchannels=4;
            %initval.chan_ID=[4 0 2 3 1];  %treat them as phase msk chr ori ter
     case 0,initval.expi='VersionTest'; %classic version test 91 cells 
            initval.pth_crop=swap_path('CD_Data_in\2016_Sandro\Testdata\');                           
            initval.pth_donut=swap_path('Dropbox\CD_Data_out\2016_Sandro\2018_SJ\BatchanalysisResults\Results_VersionTest_Rfp_2Branch_NoFlip\');              
            initval.pth_repli=swap_path('Dropbox\CD_Data_out\2016_Sandro\2018_SJ\BatchanalysisResults_Repli\VersionTest\');
            initval.pth_excelfile=swap_path('Dropbox\CD_Data_out\2016_Sandro\2018_SJ\Overview_testruns_replidata_JK_merged22.xlsx');              
            initval.chan_ID=[4 0 2 3 1];  %'OLD' treat as phase msk chr ori ter
    case 1, initval.expi='Raman1';   %5-channel Mukbef with Raman
%             initval.pth_crop=swap_path('CD_Data_in\2016_Sandro\Raman1_2019.02.04 BN2820+A22+DAPI001\'); 
    case 3, initval.expi='Test';   %5-channel Mukbef with Raman
%           initval.pth_crop=swap_path('CD_Data_in\2016_Sandro\\2820\2018-06-14-2820-DAPI\2018-06-14-2820-DAPI\');              
    case 4, initval.expi='Test20200610_MukB';   %5-channel Mukbef, hardwired paths
            initval.pth_crop='C:\Users\jkerssemakers\CD_Data_in\2016_Sandro\MukB-testdata-2019.02.04 BN2820+A22+DAPI002\';
            initval.pth_repli='C:\Users\jkerssemakers\Dropbox\CD_Data_out\2016_Sandro\2018_SJ\BatchResults_Repli\Results_Test20200610_MukB_\';
            initval.pth_excelfile='C:\Users\jkerssemakers\Dropbox\CD_Data_out\2016_Sandro\2018_SJ\Overview_testruns_replidata_JK.xlsx';             
   case 5,  initval.expi='2019.02.04_BN2820+A22+DAPI002';  
            %test for MukB symmetry, same input data as 24, remote paths, add donut data            
            initval.pth_crop=adjust_remotes('X:\tnw\BN\Shared\Raman\MukB_MatP data\Datasets_Raman\2820 (Wildtype, MukB-YFP)\2019.02.04 BN2820+A22+DAPI002\');
            initval.pth_donut=adjust_remotes('X:\tnw\BN\CD\Shared\Raman van Wee L\BatchanalysisResults\Results_BN2820_02_04B_Rfp_2Branch_Flipped\');              
            initval.pth_repli='O:\CD_Data_out\2016_Sandro\2018_SJ\BatchResults_Repli\2019.02.04_BN2820+A22+DAPI002\';
            initval.pth_excelfile='O:\CD_Data_out\2016_Sandro\2018_SJ\Overview_testruns_replidata_JK_merged22.xlsx';
            initval.channelorder=2;
            initval.numberofchannels=5;
   
%% from 23-34, run indices identical to Raman's config, transferred by JK
%
%   BN2820-Wildtype: 
%         23: 2019.01.17_BN2820+A22+42C+DAPI
%         24: 2019.02.04_BN2820+A22+DAPI001 (good images)
%         25: 2019.02.04_BN2820+A22+DAPI002
%         26: 2019.02.27_BN2820+A22+DAPI
%         30: 2019.02.27_BN2820+A22+DAPI_Tresholded_for_DNA  [1]
%         31: 2019.02.27_BN2820+A22+DAPI_Tresholded_for_foci [1]
%   BN2822-delta MatP: 
%         27: 2019.01.17_BN2822+A22+42C+DAPI.001
%         28: 2019.02.04_BN2822+A22+DAPI002
%         29: 2019.02.27_BN2822+A22+DAPI
%         32: 2019.02.27_BN2822+A22+DAPI_Tresholded_for_DNA [1]
%         33: 2019.02.27_BN2822+A22+DAPI_Tresholded_for_foci [1]
%   MukEQ strain: 
%         34: 2018_10_08_MukEQ
%
%Remarks [JK]:
%[1] [nonfunctional *] old version, only replicode output change. current version works with
% channel-specific tresholds, so these are covered by run 26 or 29
%[2] see also notes in \2021_03_02 Revisit MukB replicode 

   case 23,   initval.expi='2019.01.17_BN2820+A22+42C+DAPI';   %5-channel Mukbef with Raman
              %epth='\2820 (Wildtype, MukB-YFP)\2019.01.17 BN2820+A22+42C+DAPI\BN2820+A22+42C+DAPI.001\';              
              datapth=[drivepath 'MukB-MatP paper\MukB raw data\BN2820-Wildtype\'];
              initval.pth_crop=[datapth 'Crop Code output\BN2820+42C+DAPI\'];
              initval.pth_donut=[datapth 'Donut Code output\Results_BN2820+42C+DAPI_Rfp_SingleBranch_NoFlip\'];              
              initval.pth_repli=[datapth 'Repli Code output\BN2820+42C+DAPI\'];
              initval.pth_excelfile=[drivepath 'MukB-MatP paper\Overview_testruns_replidata.xlsx'];
   case 24,    initval.expi='2019.02.04_BN2820+A22+DAPI001';   %5-channel Mukbef with Raman
               %epth='\2820 (Wildtype, MukB-YFP)\2019.02.04 BN2820+A22+DAPI001 (good images)\'; 
              datapth=[drivepath 'MukB-MatP paper\MukB raw data\BN2820-Wildtype\'];
              initval.pth_crop=[datapth 'Crop Code output\BN2820+A22+DAPI001 (good images)\'];
              initval.pth_donut=[datapth 'Donut Code output\Results_BN2820+DAPI001_Rfp_SingleBranch_NoFlip\'];              
              initval.pth_repli=[datapth 'Repli Code output\BN2820+A22+DAPI001 (good images)\'];
              initval.pth_excelfile=[drivepath 'MukB-MatP paper\Overview_testruns_replidata.xlsx'];                
   case 25    %ori path Raman_25:
              %epth='\2820 (Wildtype, MukB-YFP)\2019.02.04 BN2820+A22+DAPI002\' 
              datapth=[drivepath 'MukB-MatP paper\MukB raw data\BN2820-Wildtype\'];
              initval.expi='2019.02.04_BN2820+A22+DAPI002'; 
              initval.pth_crop=[datapth 'Crop Code output\BN2820+A22+DAPI002\'];
              initval.pth_donut=[datapth 'Donut Code output\Results_BN2820+DAPI002_Rfp_SingleBranch_NoFlip\'];              
              initval.pth_repli=[datapth 'Repli Code output\BN2820+A22+DAPI002\'];
              initval.pth_excelfile=[drivepath 'MukB-MatP paper\Overview_testruns_replidata.xlsx'];
    case 26,  initval.expi='2019.02.27_BN2820+A22+DAPI';   
              %epth='\2820 (Wildtype, MukB-YFP)\2019.02.27 BN2820+A22+DAPI\2019.02.27 BN2820+A22+DAPI\';
              datapth=[drivepath 'MukB-MatP paper\MukB raw data\BN2820-Wildtype\'];
              initval.pth_crop=[datapth 'Crop Code output\BN2820+A22+DAPI\'];
              initval.pth_donut=[datapth 'Donut Code output\Results_BN2820+DAPI_Rfp_SingleBranch_NoFlip\'];              
              initval.pth_repli=[datapth 'Repli Code output\BN2820+A22+DAPI\'];
              initval.pth_excelfile=[drivepath 'MukB-MatP paper\Overview_testruns_replidata.xlsx'];
    case 27,  initval.expi='2019.01.17_BN2822+A22+42C+DAPI.001';   
              %epth='\2822 (Delta MatP, MukB-YFP)\2019.01.17\BN2822+A22+42C+DAPI.001\';
              datapth=[drivepath 'MukB-MatP paper\MukB raw data\BN2822-delta MatP\'];
              initval.pth_crop=[datapth 'Crop Code output\BN2822+A22+DAPI\'];
              initval.pth_donut=[datapth 'Donut Code output\Results_A22+DAPI_Rfp_SingleBranch_NoFlip\'];              
              initval.pth_repli=[datapth 'Repli Code output\BN2822+A22+DAPI\'];
              initval.pth_excelfile=[drivepath 'MukB-MatP paper\Overview_testruns_replidata.xlsx'];
    case 28 %ori path Raman_28:            
              %epth='\2822 (Delta MatP, MukB-YFP)\2019.02.04 BN2822+A22+DAPI002\';              
              datapth=[drivepath 'MukB-MatP paper\MukB raw data\BN2822-delta MatP\'];
              initval.expi='2019.02.04_BN2822+A22+DAPI002'; 
              initval.pth_crop=[datapth 'Crop code output\BN2822+A22+DAPI002\'];
              initval.pth_donut=[datapth 'Donut Code output\Results_A22+DAPI002_Rfp_SingleBranch_NoFlip\'];              
              initval.pth_repli=[datapth 'Repli Code output\BN2822+A22+DAPI002\'];
              initval.pth_excelfile=[drivepath 'MukB-MatP paper\Overview_testruns_replidata.xlsx'];
    case 29,  initval.expi='2019.02.27_BN2822+A22+DAPI';   %5-channel Mukbef with Raman
              %epth='\2822 (Delta MatP, MukB-YFP)\2019.02.27 BN2822+A22+DAPI\2019.02.27 BN2822+A22+DAPI\';
              datapth=[drivepath 'MukB-MatP paper\MukB raw data\BN2822-delta MatP\'];
              initval.pth_crop=[datapth 'Crop Code output\BN2822\'];
              initval.pth_donut=[datapth 'Donut Code output\Results_BN2822_Rfp_SingleBranch_NoFlip\'];              
              initval.pth_repli=[datapth 'Repli Code output\BN2822\'];
              initval.pth_excelfile=[drivepath 'MukB-MatP paper\Overview_testruns_replidata.xlsx'];
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
              datapth=[drivepath 'MukB-MatP paper\MukB raw data\MukEQ strain\'];
              initval.pth_crop=[datapth 'Crop Code output\2018_10_08_MukEQ\'];
              initval.pth_donut=[datapth 'Donut Code output\Results_MukEQ_Rfp_SingleBranch_NoFlip\'];              
              initval.pth_repli=[datapth 'Repli Code output\2018_10_08_MukEQ\'];
              initval.pth_excelfile=[drivepath 'MukB-MatP paper\Overview_testruns_replidata.xlsx'];           
end

       
if ismac, initval.DirSep='/';else initval.DirSep='\';end;
initval.codepth=[pwd '\'];
addpath(genpath(initval.codepth));

initval.maindatapath=strcat(initval.pth_crop, 'X050_densitydata\');
initval=get_exclusionlist(initval,batchrunindex);
if ~isdir(initval.pth_repli), mkdir(initval.pth_repli); end
save(strcat(initval.pth_repli,initval.DirSep,initval.expi, '_settings_last_run.mat'), 'initval');



function initval=get_exclusionlist(initval,batchrunindex);
%% exclusion list
switch batchrunindex
    case 0,    %'VersionTest';          epth='Testdata\';
        initval.exclusionlist=[{'100170'},{[]},{[]},{[]},{[]},{[]}];
    case 1,   %'Raman1, '5-channel Mukbef with Raman       
        initval.exclusionlist=[{'cell_002xy3'},{'cell_002xy2'},{[]},{[]},{[]},{[]}];  
      %copied from Raman (see \2021_03_02 Revisit MukB replicode notes)
    case 23,   %'Raman1, '5-channel Mukbef with Raman       
        initval.exclusionlist=[{'cell_007xy7'},{'cell_008xy1'},{'cell_008xy4'},{'cell_009xy1'},{'cell_009xy8'},{'cell_010xy5'},{'cell_010xy8'},{'cell_011xy3'},{'cell_011xy6'},{'cell_011xy5'},{'cell_011xy8'},{'cell_012xy8'},{'cell_012xy9'},{'cell_013xy2'},{'cell_013xy6'},{'cell_013xy9'},{'cell_014xy3'},{'cell_014xy5'},{'cell_015xy1'},{'cell_016xy1'},{'cell_016xy2'},{'cell_016xy3'},{'cell_016xy4'},{'cell_016xy5'},{'cell_016xy6'},{'cell_016xy7'},{'cell_017xy2'},{'cell_017xy7'},{'cell_017xy8'},{'cell_017xy9'},{'cell_018xy1'},{'cell_018xy2'},{'cell_018xy4'},{'cell_018xy6'},{'cell_018xy9'},{'cell_019xy1'},{'cell_019xy6'},{'cell_020xy4'},{'cell_020xy6'},{'cell_021xy1'},{'cell_021xy4'},{'cell_021xy7'},{'cell_022xy1'},{'cell_022xy2'},{'cell_022xy6'},{'cell_023xy1'},{'cell_023xy8'},{'cell_023xy9'},{'cell_024xy2'},{'cell_024xy1'},{'cell_025xy6'},{'cell_026xy1'},{'cell_028xy7'},{'cell_028xy8'},{'cell_029xy8'},{'cell_030xy1'},{'cell_030xy5'},{'cell_030xy8'},{'cell_031xy5'},{'cell_031xy6'},{'cell_031xy7'},{'cell_032xy5'},{'cell_032xy6'},{'cell_032xy8'},{'cell_033xy7'},{'cell_034xy5'},{'cell_035xy8'},{'cell_036xy1'},{'cell_037xy5'},{'cell_037xy6'},{'cell_037xy8'},{'cell_038xy5'},{'cell_038xy6'}];
    case 25,   %'Raman1, '5-channel Mukbef with Raman       
        initval.exclusionlist=[{'cell_009xy1'},{'cell_009xy2'},{'cell_009xy6'},{'cell_009xy7'},{'cell_010xy2'},{'cell_010xy6'},{'cell_011xy1'},{'cell_011xy2'},{'cell_011xy4'},{'cell_012xy7'},{'cell_013xy1'},{'cell_013xy2'},{'cell_014xy4'},{'cell_015xy6'},{'cell_015xy7'},{'cell_016xy1'},{'cell_016xy3'},{'cell_017xy3'},{'cell_017xy4'},{'cell_017xy6'},{'cell_018xy1'},{'cell_018xy3'},{'cell_018xy6'},{'cell_019xy6'},{'cell_019xy7'},{'cell_020xy1'},{'cell_020xy6'},{'cell_021xy1'},{'cell_021xy8'},{'cell_022xy4'},{'cell_022xy9'},{'cell_022xy3'},{'cell_025xy9'},{'cell_026xy9'},{'cell_027xy8'}];
    case 27,   %'Raman1, '5-channel Mukbef with Raman       
        initval.exclusionlist=[{'cell_009xy5'},{'cell_011xy1'},{'cell_011xy4'},{'cell_014xy1'},{'cell_014xy4'},{'cell_017xy3'},{'cell_019xy3'},{'cell_021xy4'}];
    case 28,   %'Raman1, '5-channel Mukbef with Raman       
        initval.exclusionlist=[{'cell_007xy8'},{'cell_008xy1'},{'cell_008xy2'},{'cell_008xy4'},{'cell_009xy1'},{'cell_009xy2'},{'cell_009xy5'},{'cell_010xy1'},{'cell_010xy4'},{'cell_010xy7'},{'cell_011xy6'},{'cell_012xy3'},{'cell_013xy4'},{'cell_013xy6'},{'cell_013xy7'},{'cell_014xy1'},{'cell_014xy7'},{'cell_015xy2'},{'cell_015xy3'},{'cell_015xy7'},{'cell_016xy4'},{'cell_016xy7'},{'cell_016xy8'},{'cell_017xy4'},{'cell_020xy4'},{'cell_021xy4'}];
    case 34,   %'Raman1, '5-channel Mukbef with Raman       
        initval.exclusionlist=[{'cell_002xy1'},{'cell_002xy2'},{'cell_003xy3'},{'cell_004xy2'},{'cell_005xy6'},{'cell_007xy3'},{'cell_008xy1'},{'cell_009xy1'},{'cell_009xy5'},{'cell_010xy4'},{'cell_013xy6'},{'cell_005xy1'},{'cell_005xy2'},{'cell_007xy2'},{'cell_008xy4'},{'cell_013xy3'},{'cell_014xy2'},{'cell_014xy4'},{'cell_014xy5'},{'cell_016xy4'},{'cell_016xy6'},{'cell_017xy5'},{'cell_018xy1'},{'cell_020xy4'},{'cell_020xy5'},{'cell_022xy3'},{'cell_023xy1'},{'cell_024xy1'},{'cell_025xy5'},{'cell_026xy1'},{'cell_026xy2'},{'cell_027xy2'},{'cell_026xy4'},{'cell_027xy5'},{'cell_028xy2'},{'cell_028xy5'},{'cell_029xy1'},{'cell_029xy2'},{'cell_030xy4'},{'cell_030xy2'},{'cell_030xy3'},{'cell_031xy2'},{'cell_032xy1'},{'cell_032xy4'},{'cell_034xy1'},{'cell_035xy2'},{'cell_035xy3'},{'cell_037xy1'},{'cell_037xy4'},{'cell_038xy4'},{'cell_039xy1'},{'cell_039xy4'},{'cell_041xy1'},{'cell_042xy4'},{'cell_044xy4'},{'cell_007xy7'},{'cell_008xy9'},{'cell_011xy9'},{'cell_015xy8'},{'cell_018xy8'},{'cell_019xy9'},{'cell_021xy9'},{'cell_024xy7'}];
end

function drivepath=set_drivepath;
%main paths:
if exist('O:\CD_Data_out\2016_Sandro\2018_SJ\')
    drivepath='O:\CD_Data_out\2016_Sandro\2018_SJ\';
end
if exist('N:\tnw\BN\CD\Shared\Sandro\')
    drivepath='N:\tnw\BN\CD\Shared\Sandro\';
end
if exist('N:\tnw\BN\CD\Shared\Sandro\')
    drivepath='N:\tnw\BN\CD\Shared\Sandro\';
end

