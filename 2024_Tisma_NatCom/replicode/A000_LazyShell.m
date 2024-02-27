function A000_LazyShell
%JWJK_A:-------------------------------------------------------------------
%Shell 

%Description: Shell program to run sequential steps of main code for
%replication analysis. Runs all main programs in set order; former runs can be
%re-used by setting the switches 'if 1'

%input: prior 'crop code' analysis data.

%output:  new .mat database is generated and expanded per step.

%Reference: CD lab, project Sandro, written by Jacob Kers 2018-20
%:JWJK_A-------------------------------------------------------------------

usr='Jacob';
batchrunindices=[0];               %'VersionTest'; 4-chan
% batchrunindices=[1];  %'Raman1';   %5-channel Mukbef with Raman
% batchrunindices=[23]; %initval.expi='Test'; 
batchrunindices=[24.1]; %as 24, test for MukB symmetry,  %same input data as 24, remote paths, add dont data 
batchrunindices=[23 24 25 26 27 28 29 34]; %MukB paper 2022 re-assess data
batchrunindices=[23 24 25 26 27    29 34]; %MukB paper 2022 re-assess data
batchrunindices=[23 24 25 26 27 29 34]; %MukB paper 2022 re-assess data
%batchrunindices=[25 26 27 29 34]; %MukB paper 2022 re-assess data
%batchrunindices=[27 29]; 

batchrunindices=[-2 -3];  %' repotest: -1, VersionTest: 1; 91 cells 4-chan

for ii=1:length(batchrunindices)    
    batchrunindex=batchrunindices(ii);
    disp(['batchrundindex:', num2str(batchrunindex)]);
    initval=A000_Repli_Init(batchrunindex,usr);
    
    %% basic analysis steps for replication series (in order of use)
    if 1, A001_Erase_and_Build_ResultDirs(initval);end %may require multiple tries   
    if 1, A010_WF_PerCell_AnalyzeCellShapeStandAlone(initval);end
    if 1, A013_WF_PerCell_AnalyzeSpots(initval);end
    if 1, A055_WF_2DClusterAnalysis_Standalone(initval);end
    if 0, A060_Colocalization(initval);end
    if 0, A062_MukB_Dna_overlap(initval);end
            
    %% sorting and processing steps
    %not recently checked...
    %replication series
    if 0, A005_Build_MovieList(batchrunindex,initval);end  %only use when movies
    if 0, P100_AnalyzeSpotGeometries(batchrunindex,initval);end
    if 0, P110_cluster_replisome_I_vs_D(batchrunindex,initval);end
    
end
