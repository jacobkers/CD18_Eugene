function initval=A000__WF_Get_JacobPathsandExperiments(batchrunindex)
%JWJK_B:%------------------------------------------------------------------
%Settings per experiment and initialization
%Summary: This function contains general and experiment-specific names and paths. 
%In addition, experiment-specific stettings should be originally be  loaded from an 
%'Experiment overview' excel file. To add an experiment, add entries to
%this function and the excel file; use earlier entires as format example.
%More recent experiments are run from matlab-listed settings in 'A0001_WF_ConfigPerExperiment'
%:JWJK_B-------------------------------------------------------------------
if nargin<1,batchrunindex=-1;end
if ismac, initval.DirSep='/';else initval.DirSep='\';end;

%% common settings
%runtime
initval.shortset=20E6;

%% %This setting determines from which file you get the config data. 
%historically, we use Excel, the Matlab file stores same settings
if batchrunindex>1 
    readconfig='Excel';
else
    readconfig='Matlab'; 
end


% paths
initval.codepth=pwd;
addpath(initval.codepth);
addpath(strcat(initval.codepth,initval.DirSep,'tools_donut',initval.DirSep));
cd ..;
addpath(genpath(strcat(pwd,initval.DirSep,'common_tools',initval.DirSep))); 
cd(initval.codepth);
initval.mainpath=initval.codepth; 



%% experiment labels
%code test experiments:
switch batchrunindex
    case -104.1, initval.expi='20231123_BSG5522_DAPI_2';
    case -103.3, initval.expi='20230515_BSG4610';
    case -103.2, initval.expi='20230518_BSG4610';
    case -103.1 ,initval.expi='20230129_BSG219-BSG217_comparison';
    case -102.2 ,initval.expi='Tisma_20221009_4623_DAPI';
    case -102.1 ,initval.expi='Tisma_20221009_4623';
    case -101.2 ,initval.expi='Tisma_20220803_4595_SyG';
    case -101.1 ,initval.expi='Tisma_20220802_4595_SyG';
    case -100.3 ,initval.expi='Tisma_20220614_4595';  
    case -100.2 ,initval.expi='Tisma_4595_IPTG_2h_SyG_test';   
    case -100.1 ,initval.expi='Tisma__4595_IPTG_2h_test';   
    case 0,     initval.expi='repo_test_2'; 
    case -1,    initval.expi='repo_test_N';   
    case  1,    initval.expi='VersionTest';               
end

%% get detailed info for experiment
switch readconfig
    case 'Excel';   %per Excel Table
        initval=fetch_entries_by_excel(initval);        
    case 'Matlab'%per regular matlab function 
         ThisExp=A0001_WF_ConfigPerExperiment(initval);
         initval=fetch_entries_by_matlab(ThisExp,initval);                          
end

%% inferred settings
%paths
switch initval.FlipOriention
     case 0,alignlabel=strcat(initval.alignmodus,'_NoFlip');
     case 1, alignlabel=strcat(initval.alignmodus,'_Flipped');
end
initval.sourcepath=strcat(initval.mainexperimentpath, initval.epth, initval.DirSep);
initval.maindatapath=initval.sourcepath;
initval.resultpath=strcat(initval.mainoutpath,'\BatchanalysisResults\Results_',initval.expi,'_',alignlabel,initval.DirSep);

%cells
initval.Cell_Labels=Select_data(initval);

disp(strcat('Running experiment',initval.expi)); 
disp(strcat('Config from:',readconfig));

  function initval=fetch_entries_by_matlab(ThisExp,initval);
         initval.chan_ID=ThisExp.chan_ID;        
         initval.straintype=ThisExp.straintype;
         initval.searchlabel=ThisExp.searchlabel;
         initval.alignmodus=ThisExp.alignmodus;
         initval.alignmodus_A55=ThisExp.alignmodus_A55;
         initval.cellmaskname=ThisExp.MaskName;       
         initval.FlipOriention=ThisExp.FlipCells;
         initval.PassAllCells=ThisExp.PassAllCells;  
         initval.genomelength=ThisExp.genomelength;
         percentagescale=100/initval.genomelength;
         initval.StartLabelpos=(percentagescale*ThisExp.pos_c3);
         initval.StopLabelpos=(percentagescale*ThisExp.pos_c4);
         initval.globalgappos=(percentagescale*ThisExp.pos_maingap);
         initval.Oripos=round(percentagescale*ThisExp.pos_ori); 
         initval.Difpos=round(percentagescale*ThisExp.pos_dif);
         
         %Screen parameters
         initval.Screen.CellsizeMinRadius=ThisExp.ScreenCellsizeMinRadius;
         initval.Screen.CellsizeMaxStd=ThisExp.ScreenCellsizeMaxStd;
         initval.Screen.ChromosomeSizeMinRadius=ThisExp.ScreenChromosomeSizeMinRadius;  %double for Sim
         initval.Screen.ChromosomeSizeMaxStd=ThisExp.ScreenChromosomeSizeMaxStd;     %double for Sim
         initval.Screen.ChromosomeDonutHoleDepthMin=ThisExp.ScreenChromosomeDonutHoleDepthMin;  %double for Sim
         initval.Screen.MinimalLabelAngle=ThisExp.ScreenMinimalLabelAngle; %in degs
         initval.Screen.TerOriMinDistFromCenter=ThisExp.ScreenTerOriMinDistFromCenter;  %relative from av chrom radius
         initval.Screen.CellCircularity=ThisExp.ScreenCellCircularity;
         initval.Screen.ChroCircularity=ThisExp.ScreenChroCircularity;
         initval.Screen.OriSingleSpotFraction=ThisExp.ScreenOriSingleSpotFraction;
         
         %Extras
         initval.Psf_est=ThisExp.Pointspread;
         initval.YFPLeakageCorrect=ThisExp.YFPLeakageCorrect;
         initval.SimScaleCorrect=ThisExp.SimScaleCorrect;
         initval.nmperpixel=ThisExp.NmPerPixel; 
         initval.svg_exports=ThisExp.svg_exports; 
         initval.save_svg=ThisExp.save_svg;
         initval.channelnames=ThisExp.channelnames;
         
         %analysis:--------------------------
         initval.UseMeasuredPSFforClusterAnalysis=ThisExp.UseMeasuredPSFforClusterAnalysis;
         initval.Padcurves=ThisExp.Padcurves;
         initval.spaghettisortstyle=ThisExp.spaghettisortstyle; 
         %diagnosis
         initval.AddFlippMarkers=ThisExp.AddFlippMarkers;
         
         %optional
         if ~isfield(initval, 'epth'), initval.epth=ThisExp.epth; end
         if ~isfield(initval, 'mainexperimentpath'), initval.mainexperimentpath=ThisExp.mainexperimentpath; end
         if ~isfield(initval, 'mainoutpath'), initval.mainoutpath=ThisExp.mainoutpath; end
         if ~isfield(initval, 'excelpath'), initval.excelpath=ThisExp.excelpath; end
         if ~isfield(initval, 'FlipModus'), initval.FlipModus='UseGlobalMin'; end

         
    function initval=fetch_entries_by_excel(initval);
         %%classic:
         Experiment_default.channelkey=[{'mask'}, {'ter'},  {'chromosome'}, {'ori'}, {'cell'}];              
         [numdat,textdat]= xlsread(initval.excelpath);
         Headers=textdat(1,:);
         col_Exps=find(strcmp(Headers,'ExpLabel')); %identify type column
         Exps=textdat(2:end,col_Exps);
         Exp_Row=find(strcmp(Exps,initval.expi)); 
         initval.straintype=char(textdat(Exp_Row+1,(find(strcmp(Headers,'straintype')))));
         initval.searchlabel=char(textdat(Exp_Row+1,(find(strcmp(Headers,'searchlabel')))));
         initval.alignmodus=char(textdat(Exp_Row+1,(find(strcmp(Headers,'AlignModus')))));
         initval.cellmaskname=char(textdat(Exp_Row+1,(find(strcmp(Headers,'MaskName')))));       
         initval.FlipOriention=(numdat(Exp_Row,(find(strcmp(Headers,'FlipCells')))-1));
         initval.PassAllCells=(numdat(Exp_Row,(find(strcmp(Headers,'PassAllCells')))-1));;  
         initval.genomelength=(numdat(Exp_Row,(find(strcmp(Headers,'genomelength')))-1));
         percentagescale=100/initval.genomelength;
         initval.StartLabelpos=percentagescale*(numdat(Exp_Row,(find(strcmp(Headers,'pos_c3')))-1));
         initval.StopLabelpos=percentagescale*(numdat(Exp_Row,(find(strcmp(Headers,'pos_c4')))-1));
         initval.globalgappos=percentagescale*(numdat(Exp_Row,(find(strcmp(Headers,'pos_maingap')))-1));;
         initval.Oripos=round(percentagescale*(numdat(Exp_Row,(find(strcmp(Headers,'pos_ori')))-1))); 
         initval.Difpos=round(percentagescale*(numdat(Exp_Row,(find(strcmp(Headers,'pos_dif')))-1)));

         %Screen parameters
         initval.Screen.CellsizeMinRadius=...
                (numdat(Exp_Row,(find(strcmp(Headers,'ScreenCellsizeMinRadius')))-1));
         initval.Screen.CellsizeMaxStd=...
                (numdat(Exp_Row,(find(strcmp(Headers,'ScreenCellsizeMaxStd')))-1));
         initval.Screen.ChromosomeSizeMinRadius=...
             (numdat(Exp_Row,(find(strcmp(Headers,'ScreenChromosomeSizeMinRadius')))-1));  %double for Sim
         initval.Screen.ChromosomeSizeMaxStd=...
             (numdat(Exp_Row,(find(strcmp(Headers,'ScreenChromosomeSizeMaxStd')))-1));     %double for Sim
         initval.Screen.ChromosomeDonutHoleDepthMin=...
             (numdat(Exp_Row,(find(strcmp(Headers,'ScreenChromosomeDonutHoleDepthMin')))-1));  %double for Sim
         initval.Screen.MinimalLabelAngle=...
             (numdat(Exp_Row,(find(strcmp(Headers,'ScreenMinimalLabelAngle')))-1)); %in degs
         initval.Screen.TerOriMinDistFromCenter=...
             (numdat(Exp_Row,(find(strcmp(Headers,'ScreenTerOriMinDistFromCenter')))-1));  %relative from av chrom radius
         initval.Psf_est=(numdat(Exp_Row,(find(strcmp(Headers,'Pointspread')))-1));
         initval.YFPLeakageCorrect=(numdat(Exp_Row,(find(strcmp(Headers,'YFPLeakageCorrect')))-1));
         initval.SimScaleCorrect=(numdat(Exp_Row,(find(strcmp(Headers,'SimScaleCorrect')))-1));
         initval.nmperpixel=(numdat(Exp_Row,(find(strcmp(Headers,'NmPerPixel')))-1)); 
