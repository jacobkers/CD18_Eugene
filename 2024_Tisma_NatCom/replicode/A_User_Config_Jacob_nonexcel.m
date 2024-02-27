function Experiment=A0001_WF_ConfigPerExperiment(initval);
%Experiment-specific settings to be called by 'init' file
   %---------------------------------------------------------------JWJK_C:
    %Alignment styles
    %     'Rfp_2Branch';  %alignment of branches starts at rfp label      
    %     'Rfp_SingleBranch';  %alignment of whole circle starts at rfp label
    %     'Cfp_2Branch';  %alignment of bracnhes starts at cfp label
    %     'Cfp_SingleBranch';  %alignment of whole circle starts at cfp label
    %---------------------------------------------------------------:JWJK_C
ThisExpLabel=initval.expi;

%% default settings: -------------------------------------------------- 
%name and paths
Experiment_default.Batchrunindex=0;
Experiment_default.ExpLabel='default_values';
Experiment_default.epth='testdata_in\X050_densitydata';
Experiment_default.mainexperimentpath=swap_path('CD_Data_in\2016_Sandro\');
Experiment_default.mainoutpath=swap_path('Dropbox\CD_Data_out\2016_Sandro\2018_SJ\')
Experiment_default.excelpath=strcat(Experiment_default.mainoutpath,'Overview_testruns_donutdata_JK.xlsx');
%channel ID and cell handling: 
Experiment_default.chan_ID=[2 1 3 4 5];  %treat as phase msk chr ori ter
Experiment_default.alignmodus='Rfp_2Branch';               Experiment_default.FlipCells=0;                            
%genomic parameters:--------------------------------------------------
Experiment_default.genomelength=4639;                      Experiment_default.pos_rfp=3908;                           
Experiment_default.pos_cfp=1644;                           Experiment_default.pos_maingap=1644;                       
Experiment_default.pos_ori=3924;                           Experiment_default.pos_dif=1589;
%Screening:-------------------------------------------------------
Experiment_default.ScreenCellsizeMinRadius=15;             Experiment_default.ScreenCellsizeMaxStd=8;	
Experiment_default.ScreenChromosomeSizeMinRadius=5;        Experiment_default.ScreenChromosomeSizeMaxStd=3.5;
Experiment_default.ScreenChromosomeDonutHoleDepthMin=0.3;  Experiment_default.ScreenMinimalLabelAngle=30;	
Experiment_default.ScreenCellCircularity=1.2;              Experiment_default.ScreenChroCircularity=1.2;
Experiment_default.ScreenTerOriMinDistFromCenter=0.3;      Experiment_default.PassAllCells=0;
%other:--------------------------------------------------------------
Experiment_default.straintype='type 1';                    Experiment_default.Pointspread=2.7; 
Experiment_default.SimScaleCorrect=1;                      Experiment_default.NmPerPixel=64;
Experiment_default.YFPLeakageCorrect=1;                    Experiment_default.MaskName='cellmask'; 
Experiment_default.searchlabel='c';
%% update only values that are different from default:
Experiment=Experiment_default;
switch initval.expi
    case 'VersionTest'
        %Version test classic:
        Experiment=Experiment_default;
        Experiment.Batchrunindex=1;
        Experiment.ExpLabel='VersionTest';        
        Experiment.mainexperimentpath=swap_path('CD_Data_in\2016_Sandro\');
        Experiment.epth='Testdata\X050_densitydata';
        Experiment.mainoutpath=swap_path('Dropbox\CD_Data_out\2016_Sandro\2018_SJ\')        
        Experiment.excelpath=strcat(Experiment.mainoutpath,'Overview_testruns_donutdata_JK.xlsx');
        Experiment.chan_ID=[5 1 3 4 2];  %treat as phase msk chr ori ter
    case 'repo_test'
        %Version test repo:       
        Experiment.Batchrunindex=-1;
        Experiment.ExpLabel='repo_test';
        Experiment.epth='testdata_in\X050_densitydata';
        cd .. ; 
        Experiment.mainexperimentpath=[pwd, '\'];
        Experiment.mainoutpath=[Experiment.mainexperimentpath, '\testdata_out\'];
        Experiment.excelpath=strcat(Experiment.mainoutpath,'Overview_testruns_donutdata_JK.xlsx');
        cd(initval.codepth); 
        Experiment.PassAllCells=1; 
    case 'Tisma__4595_IPTG_2h_test'
        Experiment=Experiment_default;
        Experiment.ExpLabel='Tisma__4595_IPTG_2h_test';
        Experiment.Batchrunindex=-2;
        Experiment.alignmodus='Rfp_SingleBranch';
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.7;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        if isfolder('M:\tnw\bn\cd\Shared\Tisma'), Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\'; end
        if isfolder('O:\CD_Data_in\2022_Tisma\'), Experiment.mainexperimentpath='O:\CD_Data_in\2022_Tisma\'; end
        if isfolder('D:\jkerssemakers\CD_Data_in\2022_Tisma'), Experiment.mainexperimentpath='D:\jkerssemakers\CD_Data_in\2022_Tisma\'; end
        switch 1
            case 1, Experiment.mainoutpath=swap_path('Dropbox\CD_Data_out\2022_Tisma\');
            case 2, Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        end
        Experiment.epth='20220421_4595_tests_analysis\tiffcz_4595_IPTG_2h\X050_densitydata\';
    case 'Tisma_4595_IPTG_2h_SyG_test'
        Experiment=Experiment_default;
        Experiment.ExpLabel='Tisma_4595_IPTG_2h_SyG_test';
        Experiment.Batchrunindex=-3;        
        Experiment.alignmodus='Rfp_SingleBranch';
        %Screening:-------------------------------------------------------
        Experiment.ScreenCellsizeMinRadius=12;             Experiment.ScreenCellsizeMaxStd=5;	
        Experiment.ScreenChromosomeSizeMinRadius=6;        Experiment.ScreenChromosomeSizeMaxStd=2;
        Experiment.ScreenChromosomeDonutHoleDepthMin=0.7;  Experiment.ScreenMinimalLabelAngle=-1000;	
        Experiment.ScreenCellCircularity=1.2;              Experiment.ScreenChroCircularity=1.2;
        Experiment.ScreenTerOriMinDistFromCenter=0.3;      Experiment.PassAllCells=0;
        %------------------------------------------------------------------
        if isfolder('M:\tnw\bn\cd\Shared\Tisma'), Experiment.mainexperimentpath='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\'; end
        if isfolder('O:\CD_Data_in\2022_Tisma\'), Experiment.mainexperimentpath='O:\CD_Data_in\2022_Tisma\'; end
        if isfolder('D:\jkerssemakers\CD_Data_in\2022_Tisma'), Experiment.mainexperimentpath='D:\jkerssemakers\CD_Data_in\2022_Tisma\'; end
        switch 1
            case 1, Experiment.mainoutpath=swap_path('Dropbox\CD_Data_out\2022_Tisma\');
            case 2, Experiment.mainoutpath='M:\tnw\bn\cd\Shared\Jacob\TESTdata_out\2022_Tisma\';
        end
        Experiment.epth='20220421_4595_tests_analysis\tiffcz_4595_IPTG_2h_SyG\X050_densitydata\';
end
