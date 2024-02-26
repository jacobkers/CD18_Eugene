function pathinfo=X_User_config_Jacob(expno)
%% Contents
%----------------------------------
%contains settings for the following experiments:
% test files Tisma (older runs) : 100.1, 100.2, 100.3
%   * 100.1, 4595_IPTG_2h_test_001 
%   * 100.2, 4595_IPTG_2h_SyG_test_001
%   * 100.3 20220614_4595  %suffers from bugs
% extended data (with subdirs) Tisma
%   * 101.1 : 20220802_4595_SyG: c2/3 DNA/ParB
%   * 101.2 : 20220803_4595_SyG: c2/3 DNA/ParB
% pre-submission runs MS22:
%   * 102.1 %'Tisma_20221009_4623': c2/3/4 DNA/SMC/ParB
%   * 102.2 %'Tisma_20221009_4623_DAPI' c2/3/4 DNA/SMC/ParB donut incomplete?
% rebuttal runs MS22:
%   * 103.1 %'20230129_BSG219-BSG217_comparison': DNA only
%   * 103.2 %'20230518_BSG4610': c2/3 DNA/ParB
%   * 103.3 %'20230515_BSG4610': c2/3 DNA/ParB
%rebuttal-II run
%   * 104.1 %'20231123_BSG5522_DAPI_2': c1_BF;c2_chr(DAPI) c3_ori c4_ter(bad) 

%note: '+' sign for 'donut code' '-' for crop code
%------------------------------------------------

%% general settings:
%runtime:
pathinfo.maxcellsperframe=20E9;  %limited file run for test purposes; set to 1E9 to run all data
pathinfo.newchannelnaming=1; 
pathinfo.autorun=1;
%paths and file loading options:
pathinfo.dircode = pwd; % Directory of the codes
pathinfo.dirdipimage = 'C:\Program Files\DIPimage 2.9\common\dipimage'; % Directory of dipimage    
addpath(genpath(pathinfo.dircode));
addpath(pathinfo.dirdipimage);
cd ..; addpath(genpath([pwd,'\common_tools\'])); 
pathinfo.mainpath= 'O:\CD_Data_in\2016_Sandro\';

%how to pick the images:
%either by a specified plane index, or by for example max-projection

%image handling:
pathinfo.cropedge=5;        % cropping space in all directions
pathinfo.mincellsize=500;    % in image pixels 
pathinfo.maxcellsize=10000; 
pathinfo.numberofcolours=4; %historical phase chro ori ter


%% specific settings
        if nargin<1, expno=1;end
          switch expno %
              case 104.1 %'20231123_BSG5522_DAPI_2'[RebuttalII)
                    pathinfo.numberofcolours=4; %%phase chro ori ter
                    pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                    pathinfo.experimentpath=[pathinfo.mainpath '\20231123_BSG5522_DAPI_2\'];
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    pathinfo.namebases=[{'BSG5522_DAPI_2h45min_001'}, ...
                                        {'BSG5522_DAPI_2h45min_002'}, ...
                                        {'BSG5522_DAPI_2h45min_003'}, ...
                                        {'BSG5522_DAPI_2h45min_004'}, ...
                                        {'BSG5522_DAPI_2h45min_005'}, ...
                                        {'BSG5522_DAPI_2h45min_006'}, ...
                                        {'BSG5522_DAPI_2h45min_007'}, ...
                                        ];
                    N_namebases=length(pathinfo.namebases);
                    zplane_1chan=[8 8 8 8 8 8 8]';  %manual imageJ CHECKIT'          
                    pathinfo.txycz_base='xy1c1z01';        
                    pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
                      % Center plane of z-scan channels:
                    pathinfo.centerplane = repmat(zplane_1chan,1,pathinfo.numberofcolours);  %three fluo plus phase               
                    pathinfo.limitzplanes=1;
                    pathinfo.autorun=1;
                    pathinfo.mincellsize=100; 
                    pathinfo.maxcellsize=15000; 
              case 103.3 %'20230515_BSG4610'
                    pathinfo.numberofcolours=4; %%phase chro parB SMC?
                    pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                    pathinfo.experimentpath=[pathinfo.mainpath '\20230515_BSG4610\'];
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    pathinfo.namebases=[{'BSG4610_001'}, {'BSG4610_002'}, ...
                                        {'BSG4610_003'}, {'BSG4610_004'},...
                                        {'BSG4610_005'}, {'BSG4610_006'}, ...
                                        {'BSG4610_007'}, {'BSG4610_008'}, ...
                                        {'BSG4610_009'}, {'BSG4610_010'}, ...
                                        {'BSG4610_011'}, {'BSG4610_012'}, ...
                                        {'BSG4610_013'}, {'BSG4610_014'}];
                    N_namebases=length(pathinfo.namebases);
                    zplane_1chan=[6 6 6 6 6 6 5 5 5 5 5 5 5 5]';  %manual imageJ CHECKIT'          
                    pathinfo.txycz_base='xy1c1z01';        
                    pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
                      % Center plane of z-scan channels:
                    pathinfo.centerplane = repmat(zplane_1chan,1,pathinfo.numberofcolours);  %three fluo plus phase               
                    pathinfo.limitzplanes=1;
                    pathinfo.autorun=1;
                    pathinfo.mincellsize=100; 
                    pathinfo.maxcellsize=15000; 
              case 103.2 %'20230518_BSG4610'
                    pathinfo.numberofcolours=4; %%phase chro parB SMC?
                    pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                    pathinfo.experimentpath=[pathinfo.mainpath '\20230518_BSG4610\'];
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    pathinfo.namebases=[{'BSG4610_001'}, {'BSG4610_002'}, ...
                                        {'BSG4610_003'}, {'BSG4610_004'},...
                                        {'BSG4610_005'}, {'BSG4610_006'}, ...
                                        {'BSG4610_007'}, {'BSG4610_008'}, ...
                                        {'BSG4610_009'}, {'BSG4610_010'}, ...
                                        {'BSG4610_011'}, {'BSG4610_012'}, ...
                                        {'BSG4610_013'}];
                    N_namebases=length(pathinfo.namebases);
                    zplane_1chan=[6 6 6 6 7 7 7 7 7 7 7 6 6]';  %manual imageJ CHECKIT'           
                    pathinfo.txycz_base='xy1c1z01';                    
                    pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
                    % Center plane of z-scan channels:
                    pathinfo.centerplane = repmat(zplane_1chan,1,pathinfo.numberofcolours);  %three fluo plus phase               
                    pathinfo.limitzplanes=1;
                    pathinfo.autorun=1;
                    pathinfo.mincellsize=100; 
                    pathinfo.maxcellsize=15000;  
              case 103.1 %'20230129_BSG219-BSG217_comparison'
                    pathinfo.numberofcolours=3; %phase chro empty 
                    pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                    pathinfo.experimentpath=[pathinfo.mainpath '\20230129_BSG219-BSG217_comparison\'];
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    switch 2
                        case 1
                            pathinfo.namebases=[{'BSG217_Xyl60_001'}, {'BSG217_Xyl60_002'}, ...
                                                {'BSG219_noXyl_001'}, {'BSG219_noXyl_002'},...
                                                {'BSG219_noXyl_003'}, {'BSG219_noXyl_004'}, ...
                                                {'BSG219_Xyl60_001'}, {'BSG219_Xyl60_002'},...
                                                {'BSG219_Xyl60_003'}, {'BSG219_Xyl60_004'}];
                            zplane_1chan=[8 8 6 6 6 6 6 6 6 6]';  %manual imageJ'
                        case 2
                            pathinfo.namebases=[{'BSG219_Xyl60_001'}, {'BSG219_Xyl60_002'},...
                                                {'BSG219_Xyl60_003'}, {'BSG219_Xyl60_004'}];
                            zplane_1chan=[6 6 6 6]';  %manual imageJ'            
                    end      
                    pathinfo.txycz_base='xy1c1z01';
                    
                    pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
                    % Center plane of z-scan channels:
                    N_namebases=length(pathinfo.namebases);
                    
                    pathinfo.centerplane = repmat(zplane_1chan,1,pathinfo.numberofcolours);  %three fluo plus phase               
                    pathinfo.limitzplanes=1;
                    pathinfo.autorun=1;
                    pathinfo.mincellsize=100; 
                    pathinfo.maxcellsize=15000;  
                case 102.1 %'Tisma_20221009_4623'
                    pathinfo.numberofcolours=4; %phase chro parB SMC
                    pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                    pathinfo.experimentpath=[pathinfo.mainpath '\20221009_4623\'];
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    pathinfo.namebases=[{'BSG4623_001'}, {'BSG4623_002'}, {'BSG4623_003'}, {'BSG4623_004'},...
                                        {'BSG4623_005'}, {'BSG4623_006'}, {'BSG4623_007'}, {'BSG4623_008'},...
                                        {'BSG4623_009'}, {'BSG4623_010'}, {'BSG4623_011'}, {'BSG4623_012'}];                                                  
                    pathinfo.txycz_base='xy1c1z02';
                    pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
                    % Center plane of z-scan channels:
                    N_namebases=length(pathinfo.namebases);
                    zplane_1chan=[5 4 3 2 4 5 2 5 5 5 5 5]';  %by Tisma 'zplanes.txt'
                    pathinfo.centerplane = repmat(zplane_1chan,1,pathinfo.numberofcolours);  %three fluo plus phase               
                    pathinfo.limitzplanes=1;
                    pathinfo.autorun=1;
                    pathinfo.mincellsize=100; 
                    pathinfo.maxcellsize=15000; 
          
                case 102.2 %'Tisma_20221009_4623_DAPI'      
                    pathinfo.numberofcolours=4; %phase chro parB SMC  
                    pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                    pathinfo.experimentpath=[pathinfo.mainpath '\20221009_4623_DAPI\'];
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    pathinfo.namebases=[{'BSG4623_001_DAPI'}];            
                    pathinfo.txycz_base='xy1c1z02';
                    pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
                    % Center plane of z-scan channels 1-3 = [brightfield, yfp, rfp]: 
                    %set per number of namebase
                    zplane_1chan=zeros(1,1)+5;
                    pathinfo.centerplane = repmat(zplane_1chan,1,pathinfo.numberofcolours);                    
                    pathinfo.limitzplanes=1;
                    pathinfo.autorun=1;
                    pathinfo.mincellsize=100; 
                    pathinfo.maxcellsize=15000;  

               %% extended data Tisma : 101.1,101.2
                case 101.1
                    pathinfo.numberofcolours=3; %phase chro parB
                    pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                    pathinfo.experimentpath=[pathinfo.mainpath '\20220802_4595_SyG\'];
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    pathinfo.namebases=[{'4595_SyG_001'}, {'4595_SyG_002'}, {'4595_SyG_003'}, {'4595_SyG_004'}, {'4595_SyG_005'},...
                                        {'4595_SyG_006'}, {'4595_SyG_007'}, {'4595_SyG_008'}, {'4595_SyG_009'}, {'4595_SyG_010'},...
                                        {'4595_SyG_011'}, {'4595_SyG_012'}, {'4595_SyG_013'}, {'4595_SyG_014'}, {'4595_SyG_015'},...
                                        {'4595_SyG_016'}, {'4595_SyG_017'}, {'4595_SyG_018'}, {'4595_SyG_019'}, {'4595_SyG_020'},...
                                        {'4595_SyG_021'}, {'4595_SyG_022'}, {'4595_SyG_023'}, {'4595_SyG_024'}, {'4595_SyG_025'},...
                                        {'4595_SyG_026'}, {'4595_SyG_027'}];                 
                    pathinfo.txycz_base='xy1c1z02';
                    pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
                    % Center plane of z-scan channels 1-3 = [brightfield, yfp, rfp]: 
                    %set per number of namebase
                    N_namebases=length(pathinfo.namebases);
                    zplane_1chan=[ zeros(N_namebases,1)+6];
                    pathinfo.centerplane = repmat(zplane_1chan,1,pathinfo.numberofcolours);                    
                    pathinfo.limitzplanes=1;
                    pathinfo.autorun=1;
                    
                 case 101.2
                     pathinfo.numberofcolours=3; %phase chro parB
                     pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                    pathinfo.experimentpath=[pathinfo.mainpath '\20220803_4595_SyG\'];
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    pathinfo.namebases=[{'4595_SyG_001'}, {'4595_SyG_002'}, {'4595_SyG_003'}, {'4595_SyG_004'}, {'4595_SyG_005'},...
                                        {'4595_SyG_006'}, {'4595_SyG_007'}, {'4595_SyG_008'}, {'4595_SyG_009'}, {'4595_SyG_010'},...
                                        {'4595_SyG_011'}, {'4595_SyG_012'}, {'4595_SyG_013'}, {'4595_SyG_014'}, {'4595_SyG_015'},...
                                        {'4595_SyG_016'}, {'4595_SyG_017'}, {'4595_SyG_018'}, {'4595_SyG_019'}, {'4595_SyG_020'},...
                                        {'4595_SyG_021'}, {'4595_SyG_022'}, {'4595_SyG_023'}, {'4595_SyG_024'}, {'4595_SyG_025'},...
                                        {'4595_SyG_026'}, {'4595_SyG_027'}, {'4595_SyG_028'}];  
                    
                                    
                    %pathinfo.namebases=[{'4595_SyG_001'} ,{'4595_SyG_002'} ];  
                    pathinfo.txycz_base='xy1c1z02';
                    pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];                   
                    % Center plane of z-scan channels 1-3 = [brightfield, yfp, rfp]:
                    %Movie 1-12:z = 8/11;Movie 13-17:z = 5/11;Movie 18-28:z = 6/11:
                    zplane_1chan=[zeros(12,1)+8; zeros(5,1)+5 ; zeros(11,1)+6 ];
                    pathinfo.centerplane = repmat(zplane_1chan,1,3);   
                    pathinfo.limitzplanes=1;
                    pathinfo.autorun=1;
                    
                    
                    
               %%  Test files Tisma : 100.1, 100.2, 100.3:
                      % 4595_IPTG_2h_SyG_test_001xy1c2_cmle_t00 [deconvolved, donuts/ U-shapes]
                        % 4595_IPTG_2h_SyG_test_001xy1c1: phase
                        % 4595_IPTG_2h_SyG_test_001xy1c2: raw dna
                        % 4595_IPTG_2h_SyG_test_001xy1c3: ParB labeled ori?
                        % 4595_IPTG_2h_SyG_test_001xy1c4: phase #2
                 case 100.3   %Test files Tisma III 
                 switch 2                 
                     case 1, pathinfo.mainpath='V:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                     case 2, pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                     case 3, pathinfo.mainpath='O:\CD_Data_in\2022_Tisma\';
                     case 4, pathinfo.mainpath='C:\Users\jkerssemakers\CD_Data_in\2022_Tisma\';
                 end
                 pathinfo.experimentpath=[pathinfo.mainpath '\20220614_4595\'];
                 pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                 pathinfo.centerplane = [2 2 2]; % Center plane of z-scan in different channels; 1-3 = [brightfield, yfp, rfp]
                 pathinfo.limitzplanes=1;
                 pathinfo.namebases=[{'4595_SyG_001'}, {'4595_SyG_002'}, {'4595_SyG_003'}, {'4595_SyG_004'}, {'4595_SyG_005'}];                 
                 pathinfo.txycz_base='xy1c1z1';
                 pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
                 pathinfo.autorun=1;
               case 100.2   %Test files Tisma II 
                 switch 2                 
                     case 1, pathinfo.mainpath='V:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                     case 2, pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                     case 3, pathinfo.mainpath='O:\CD_Data_in\2022_Tisma\';
                 end
                 pathinfo.experimentpath=[pathinfo.mainpath '\20220421_4595_tests_analysis\tiffcz_4595_IPTG_2h_SyG\'];
                 pathinfo.centerplane = [6 6 6]; % Center plane of z-scan in different channels; 1-3 = [brightfield, yfp, rfp]
                 pathinfo.limitzplanes=1;
                 pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                 pathinfo.namebases=[{'4595_IPTG_2h_SyG_test_001'}];                 
                 pathinfo.txycz_base='xy1c1z6';
                 pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif']; 
              case 100.1   %Test files Tisma remote 
                 switch 2                 
                     case 1, pathinfo.mainpath='V:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                     case 2, pathinfo.mainpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
                     case 3, pathinfo.mainpath='O:\CD_Data_in\2022_Tisma\';
                     case 4, pathinfo.mainpath='D:\jkerssemakers\CD_Data_in\2022_Tisma\';
                 end 
                pathinfo.experimentpath=[pathinfo.mainpath '\20220421_4595_tests_analysis\tiffcz_4595_IPTG_2h\'];
                pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                pathinfo.centerplane = [6 6 6]; % Center plane of z-scan in different channels; 1-3 = [brightfield, yfp, rfp]
                pathinfo.limitzplanes=1;
                pathinfo.namebases=[{'4595_IPTG_2h_test_001'}];
                pathinfo.txycz_base='xy1c1z6';
                pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
            case -1   %  repo test many cells           
                    cd(pathinfo.dircode); cd .. ; 
                    pathinfo.experimentpath=[pwd,'\testdata_in\testNcells\']; 
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    cd(pathinfo.dircode);
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    pathinfo.namebases=[{'4595_SyG_001'}];                 
                    pathinfo.txycz_base='xy1c1z02';
                    pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
                    % Center plane of z-scan channels 1-3 = [brightfield, yfp, rfp]: 
                    %set per number of namebase
                    N_namebases=length(pathinfo.namebases);
                    zplane_1chan=[ zeros(N_namebases,1)+6];
                    pathinfo.centerplane = repmat(zplane_1chan,1,pathinfo.numberofcolours);                    
                    pathinfo.limitzplanes=1;
                    pathinfo.autorun=1;    
            case 0   %repo-based test 2 cells
                    cd(pathinfo.dircode); cd .. ; 
                    pathinfo.experimentpath=[pwd,'\testdata_in\test2cells\']; 
                    pathinfo.dirraw =[pathinfo.experimentpath, '\tiffs\']; % Raw images exported from the microscope
                    cd(pathinfo.dircode);
                    %center z plane per channel; 
                    %1-3 = [brightfield, yfp,rfp]:
                    pathinfo.numberofcolours=3; 
                    pathinfo.centerplane = [6 6 6]; %if a single row, all subdirs same 
                    pathinfo.limitzplanes=0; 
                    %format of tif file names:
                    pathinfo.namebases=[{'2179-initiations_A1-1'}];
                    pathinfo.txycz_base='t01xy1c1z1';
                    pathinfo.txycz_template=[char(pathinfo.namebases{1}) pathinfo.txycz_base '.tif'];
           
          end
             
            end       