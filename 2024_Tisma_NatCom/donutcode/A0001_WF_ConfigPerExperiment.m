function Experiment=A0001_WF_ConfigPerExperiment(initval);
%JWJK_B:%------------------------------------------------------------------
% Configuration Per Experiment
% Experiment-specific settings to be called by 'init' file. General
% framework sets default settings, which overrides per specific experiment
% if wished.
    %Alignment styles
    %     'c3_2Branch';  %alignment of branches starts at c3 label      
    %     'c3_SingleBranch';  %alignment of whole circle starts at c3 label
    %     'c4_2Branch';  %alignment of bracnhes starts at c4 label
    %     'c4_SingleBranch';  %alignment of whole circle starts at c4 label
%:JWJK_B-------------------------------------------------------------------

%----------------------------------
%contains settings for the following experiments:
% test files Tisma (older runs) : 100.1, 100.2, 100.3
%   * 100.1, 4595_IPTG_2h_test_001 
%   * 100.2, 4595_IPTG_2h_SyG_test_001
%   * 100.3 20220614_4595
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
%   * 104.1 %'20231123_BSG5522_DAPI_2'; c2/3/4 DNA/ori/ter
%note: '+' sign for 'donut code' '-' for crop code
%------------------------------------------------



ThisExpLabel=initval.expi;

%% default settings: -------------------------------------------------- 
%name and paths
Experiment_default.Batchrunindex=0;
Experiment_default.ExpLabel='default_values';
Experiment_default.epth='testdata_in\X050_densitydata';
Experiment_default.mainexperimentpath=swap_path('CD_Data_in\2016_Sandro\');
Experiment_default.mainoutpath=swap_path('Dropbox\CD_Data_out\2016_Sandro\2018_SJ\');
Experiment_default.excelpath=strcat(Experiment_default.mainoutpath,'Overview_testruns_donutdata_JK.xlsx');
%channel ID, number of colors and cell handling: 
Experiment_default.chan_ID=[1 0 2 3 4];  %treatment key, as phase msk chr label1 label2
Experiment_default.channelnames=[{'mask'}, {'phase'}, {'DNA'}, {'ori'}, {'ter'}];  
                                %in order c0 c1 c2 c3 c4
Experiment_default.alignmodus_A55='c4';  %alignment of profiles in A55
Experiment_default.alignmodus='c3_2Branch';               Experiment_default.FlipCells=0; 
Experiment_default.FlipModus='UseGlobalMin'; 
%genomic parameters:--------------------------------------------------
Experiment_default.genomelength=4639;                      Experiment_default.pos_c3=3908;                           
Experiment_default.pos_c4=1644;                           Experiment_default.pos_maingap=1644;                       
Experiment_default.pos_ori=3924;                           Experiment_default.pos_dif=1589;
%Screening:-------------------------------------------------------
Experiment_default.ScreenCellsizeMinRadius=15;             Experiment_default.ScreenCellsizeMaxStd=8;	
Experiment_default.ScreenChromosomeSizeMinRadius=5;        Experiment_default.ScreenChromosomeSizeMaxStd=3.5;
Experiment_default.ScreenChromosomeDonutHoleDepthMin=0.3;  Experiment_default.ScreenMinimalLabelAngle=0;	
Experiment_default.ScreenCellCircularity=1.2;              Experiment_default.ScreenChroCircularity=1.2;
Experiment_default.ScreenTerOriMinDistFromCenter=0.3;      Experiment_default.ScreenOriSingleSpotFraction=0.0;
Experiment_default.PassAllCells=0;
%other:--------------------------------------------------------------
Experiment_default.straintype='type 1';                    Experiment_default.Pointspread=2.7; 
Experiment_default.SimScaleCorrect=1;                      Experiment_default.NmPerPixel=64;
Experiment_default.YFPLeakageCorrect=1;                    Experiment_default.MaskName='cellmask'; 
Experiment_default.searchlabel='c';
Experiment_default.svg_exports=[{'xxx'}];
Experiment_default.save_svg=1;
%analysis:--------------------------
Experiment_default.UseMeasuredPSFforClusterAnalysis=0;
Experiment_default.Padcurves=50;  %if >0, ...% periodic elongation of density curves 
Experiment_default.spaghettisortstyle='Peak'; %'Similarity'; 'Peak' %'Std';  'Imposed' 'Minima' 'Similarity'
%diagnosis
Experiment_default.AddFlippMarkers=0;   %flip check


%% update only values that are different from default:
Experiment=Experiment_default;
switch initval.expi
    case '20231123_BSG5522_DAPI_2'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-104.1;
        Experiment.ExpLabel='20231123_BSG5522_DAPI_2';        
        Experiment.alignmodus='c3_SingleBranch';  
        Experiment.alignmodus_A55='c3';  %alignment of profiles in A55
        Experiment.FlipCells=1; 
        Experiment.chan_ID=[1 0 2 3 4];  %treat as phase msk chr  ori ter
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=4;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.65; Experiment.ScreenMinimalLabelAngle=0;	
        Experiment.ScreenCellCircularity=1.4;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.4;      Experiment.ScreenOriSingleSpotFraction=0.1; 
        Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
        Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        Experiment.epth='20231123_BSG5522_DAPI_2\X050_densitydata\';
    case '20230515_BSG4610'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-103.3;
        Experiment.ExpLabel='20230515_BSG4610';        
        Experiment.alignmodus='c3_SingleBranch';
        Experiment.FlipCells=1; 
        Experiment.chan_ID=[1 0 2 3 4];  %treat as phase msk chr  SMC parB
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.9;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.ScreenOriSingleSpotFraction=0.6; 
        Experiment.PassAllCells=1;
        %------------------------------------------------------------------
        Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
        Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        Experiment.epth='20230515_BSG4610\X050_densitydata\';    
            
            
        case '20230518_BSG4610'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-103.2;
        Experiment.ExpLabel='20230518_BSG4610';        
        Experiment.alignmodus='c3_SingleBranch';
        Experiment.FlipCells=1; 
        Experiment.chan_ID=[1 0 2 3 4];  %treat as phase msk chr  SMC parB
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;      Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.9;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.ScreenOriSingleSpotFraction=0.0; 
        Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
        Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        Experiment.epth='20230518_BSG4610\X050_densitydata\';    
            
            
        case '20230129_BSG219-BSG217_comparison'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-103.1;
        Experiment.ExpLabel='20230129_BSG219';        
        Experiment.alignmodus='c3_SingleBranch';
        Experiment.FlipCells=1; 
        Experiment.chan_ID=[1 0 2 2 2];  %treat as phase msk chr  SMC parB
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.9;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.ScreenOriSingleSpotFraction=0.0; 
        Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
        Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        Experiment.epth='20230129_BSG219-BSG217_comparison\X050_densitydata\';    
        case 'Tisma_20221009_4623_DAPI'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-102.2;
        Experiment.ExpLabel='20221009_4623';        
        Experiment.alignmodus='c3_SingleBranch';
        Experiment.FlipCells=1; 
        Experiment.chan_ID=[1 0 2 3 4];  %treat as phase msk chr  SMC parB
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.9;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.ScreenOriSingleSpotFraction=0.6; 
        Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
        Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        Experiment.epth='20221009_4623_DAPI\X050_densitydata\';
    case 'Tisma_20221009_4623'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-102.1;
        Experiment.ExpLabel='20221009_4623_2subdirs';        
        Experiment.alignmodus='c3_SingleBranch';  %=SMC....
        Experiment.FlipCells=1; 
        Experiment.chan_ID=[1 0 2 3 4];  %treat as phase msk chr  SMC parB
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=10;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=3;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.9;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.ScreenOriSingleSpotFraction=0.6; 
        Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
        Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        Experiment.epth='20221009_4623\X050_densitydata\';

    case 'Tisma_20220803_4595_SyG'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-101.2;
        Experiment.ExpLabel='20220803_4595_SyG';        
        Experiment.alignmodus='c3_SingleBranch';
        Experiment.FlipCells=1; 
        Experiment.chan_ID=[1 0 2 4 3];  %treat as phase msk chr parB parB
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.9;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.ScreenOriSingleSpotFraction=0.6; 
        Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
        Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        Experiment.epth='20220803_4595_SyG\X050_densitydata\';
        Experiment.svg_exports=[{'cell_031xy1'},{'cell_101xy1'},{'cell_132xy1'},...
                                {'cell_722xy1'},{'cell_723xy1'},{'cell_773xy1'}];
      
    case 'Tisma_20220802_4595_SyG'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-101.1;
        Experiment.ExpLabel='20220802_4595_SyG';        
        Experiment.alignmodus='c3_SingleBranch';
        Experiment.FlipCells=1;
        Experiment.chan_ID=[1 0 2 4 3];  %treat as phase msk chr parB parB
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.9;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.ScreenOriSingleSpotFraction=0.6; 
        Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\';
        Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        Experiment.epth='20220802_4595_SyG\X050_densitydata\';
        Experiment.svg_exports=[{'cell_074xy1'},{'cell_125xy1'},{'cell_505xy1'},...
                                {'cell_573xy1'},{'cell_605xy1'},{'cell_625xy1'},...
                                {'cell_627xy1'},{'cell_658xy1'},{'cell_685xy1'},...
                                {'cell_795xy1'},{'cell_854xy1'},{'cell_988xy1'},...
                                {'cell_1031xy1'},{'cell_1242xy1'},{'cell_1441xy1'},...
                                {'cell_1468xy1'},{'cell_1484xy1'}];
   
    case 'Tisma_20220614_4595'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-100.3;
        Experiment.ExpLabel='Tisma_20220614_4595';        
        Experiment.alignmodus='c3_SingleBranch';
        Experiment.FlipCells=1; 
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.9;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.ScreenOriSingleSpotFraction=0.6; 
        Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        if isfolder('M:\tnw\bn\cd\Shared\Tisma\'), Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\\'; end
        if isfolder('C:\Users\jkerssemakers\CD_Data_in\2022_Tisma'), Experiment.mainexperimentpath='C:\Users\jkerssemakers\CD_Data_in\2022_Tisma\'; end
        switch 1
            case 1, Experiment.mainoutpath=swap_path('Dropbox\CD_Data_out\2022_Tisma\');
            case 2, Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        end
        Experiment.epth='20220614_4595\X050_densitydata\';
    
     case 'Tisma_4595_IPTG_2h_SyG_test'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-100.2; 
        Experiment.ExpLabel='Tisma_4595_IPTG_2h_SyG_test';               
        Experiment.alignmodus='c3_SingleBranch';
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.7;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.ScreenOriSingleSpotFraction=0.6;
        Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        if isfolder('M:\tnw\bn\cd\Shared\Tisma'), Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\'; end
        if isfolder('O:\CD_Data_in\2022_Tisma\'), Experiment.mainexperimentpath='O:\CD_Data_in\2022_Tisma\'; end
        if isfolder('D:\jkerssemakers\CD_Data_in\2022_Tisma'), Experiment.mainexperimentpath='D:\jkerssemakers\CD_Data_in\2022_Tisma\'; end
        switch 1
            case 1, Experiment.mainoutpath=swap_path('Dropbox\CD_Data_out\2022_Tisma\');
            case 2, Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        end
        Experiment.epth='20220421_4595_tests_analysis\tiffcz_4595_IPTG_2h_SyG\X050_densitydata\';
    case 'Tisma__4595_IPTG_2h_test'
        Experiment=Experiment_default;
        Experiment.Batchrunindex=-100.1;
        Experiment.ExpLabel='Tisma__4595_IPTG_2h_test';        
        Experiment.alignmodus='c3_SingleBranch';
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.7;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.ScreenOriSingleSpotFraction=0.6;
        Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        if isfolder('M:\tnw\bn\cd\Shared\Tisma'), Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\'; end
        if isfolder('O:\CD_Data_in\2022_Tisma\'), Experiment.mainexperimentpath='O:\CD_Data_in\2022_Tisma\'; end
        if isfolder('D:\jkerssemakers\CD_Data_in\2022_Tisma'), Experiment.mainexperimentpath='D:\jkerssemakers\CD_Data_in\2022_Tisma\'; end
        switch 1
            case 1, Experiment.mainoutpath=swap_path('Dropbox\CD_Data_out\2022_Tisma\');
            case 2, Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        end
        Experiment.epth='20220421_4595_tests_analysis\tiffcz_4595_IPTG_2h\X050_densitydata\';
    case 'VersionTest'
        %Version test classic:
        Experiment=Experiment_default;
        Experiment.Batchrunindex=1;
        Experiment.ExpLabel='VersionTest';        
        Experiment.mainexperimentpath=swap_path('CD_Data_in\2016_Sandro\');
        Experiment.epth='Testdata\X050_densitydata';
        Experiment.mainoutpath=swap_path('Dropbox\CD_Data_out\2016_Sandro\2018_SJ\')        
        Experiment.excelpath=strcat(Experiment.mainoutpath,'Overview_testruns_donutdata_JK.xlsx');
        Experiment.chan_ID=[4 0 2 3 1];  %treat as phase msk chr ori ter
        Experiment.svg_exports=[{'100190'},{'100227'},{'100233'}];
    case 'repo_test_N'
        %Version test repo:       
        Experiment.Batchrunindex=-1;
        Experiment.ExpLabel='repo_test_N';        
        cd .. ; 
        Experiment.mainexperimentpath=[pwd, '\'];
        Experiment.epth='\testdata_in\testNcells\X050_densitydata';
        Experiment.mainoutpath=[Experiment.mainexperimentpath, '\testdata_out\test_out_Ncells\'];
        Experiment.excelpath=[];
        cd(initval.codepth); 
        Experiment.PassAllCells=1; 
    case 'repo_test_2'
        %Version test repo:       
        Experiment.Batchrunindex=0;
        Experiment.ExpLabel='repo_test_2';      
        cd .. ; 
        Experiment.mainexperimentpath=[pwd, '\'];  
        Experiment.epth='\testdata_in\test2cells\X050_densitydata';
        Experiment.mainoutpath=[Experiment.mainexperimentpath, '\testdata_out\test_out_2cells\'];
        Experiment.excelpath=[];
        cd(initval.codepth); 
        Experiment.PassAllCells=1;      
end
