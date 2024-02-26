function A000_WF__ClickAndGo
%JWJK_A:----[add JWJK_:ABCorC*]----------------------------------------------------
%ClickAndGo shell
%Description: Shell program for batch analysis of chromatin imaging
%experiments. Check or uncheck the following sub-sections (a first
%run of every step is necessary)
    % * A010_WF_PerCell_AnalyzeAllColors
    % * A015_WF_PerCell_Okaying(batchrunindex)
    % * A020_WF_PerCell_ChromosomeAlignment(batchrunindex)
    % * A030_WF_PlotSpaghettiCurves(batchrunindex)
    % * A040_WF_PlotGapsAndPeaks(batchrunindex)
    % * A050_WF_2DClusterAnalysis(batchrunindex)
    % * A055_WF_ReSampleContours(batchrunindex)
%Input: none
%Output: various
%References: Jacob Kerssemakers, Cees Dekker Lab, Delft
%:JWJK_A-----------------------------------------------------
%
%% choose experiment
%repo_test :0 (2 cells) and -1 (many cells)
%BatchrunExpArray=[-1 1];



%2023 rebuttal synchronization runs:
BatchrunExpArray = [ ...
    -104.1,...  %'20231123_BSG5522_DAPI_2'; 
    -103.3,...  % '20230515_BSG4610';
    -103.2,...  % '20230518_BSG4610';
    -103.1,...  % '20230129_BSG219-BSG217_comparison';
    ..., %-102.2 ,...  % 'Tisma_20221009_4623_DAPI';
    -102.1 ,...  % 'Tisma_20221009_4623';
    -101.2 ,...  % 'Tisma_20220803_4595_SyG';
    -101.1 ,...  % 'Tisma_20220802_4595_SyG';
    -100.3 ,...  % 'Tisma_20220614_4595';  
     ];
BatchrunExpArray = [ ...
    -102.1 ,...  % 'Tisma_20221009_4623';
    -101.2 ,...  % 'Tisma_20220803_4595_SyG';
    -101.1 ,...  % 'Tisma_20220802_4595_SyG';
    -100.3 ,...  % 'Tisma_20220614_4595';  
     ];
BatchrunExpArray = [-102.1 ]% 'Tisma_20221009_4623';


 %for the dryrun:
dryrun=0; %just to check if all data is ready to run
for ii=1:length(BatchrunExpArray)
    batchrunindex=BatchrunExpArray(ii);
    if 1, A010_WF_PerCell_AnalyzeAllColors(batchrunindex); end
    if 1, A015_WF_PerCell_Okaying(batchrunindex,dryrun); end        
    if 1, A020_WF_PerCell_ChromosomeAlignment(batchrunindex), end
    %not used in Tisma Project:------------------------
    if 0, A030_WF_PlotSpaghettiCurves(batchrunindex), end 
    if 0, A040_WF_PlotGapsAndPeaks(batchrunindex), end
    %---------------------------------------------
    if 1, A050_WF_2DClusterAnalysis(batchrunindex);end
    if 1, A051_WF_Work_Clusters(batchrunindex), end
    if 1, A055_WF_Work_Demographs(batchrunindex), end
    
end

