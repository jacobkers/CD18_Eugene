function X0000_AutoDataRunShell
%this shell is intended for quick running of the crop code sequence on new
%and existing data. Note that it may overwrite results of manual existing runs!


%% choose user and experimentindex
usr='Jacob';
expstotestrun=[0 -1]; %repository test on 2 cells (code check)
expstotestrun=[0]; %repository test on 2 cells (code check)  
expstotestrun=[-1]; %repository test on 141 cells (data check) 
%data runs (see config)
expstotestrun=[100.1 100.2 100.3];  %test runs; DNA/ParB, some bugs
expstotestrun=[101.1 101.2 100.3];  %extended subdirs; DNA/ParB
expstotestrun=[102.1 102.2];  %DNA/SMC/ParB
expstotestrun=[103.1 103.2 103.3];  %rebuttal runs; DNA, DNA/Parb, DNA/Parb
%all:
expstotestrun=[
    [101.1 101.2],...;  %extended subdirs; DNA/ParB
    [102.1 102.2],...;  %DNA/SMC/ParB
    [103.1 103.2 103.3]];  %rebuttal runs; DNA, DNA/Parb, DNA/Parb
expstotestrun=[
    [102.1 102.2],...;  %DNA/SMC/ParB
    [103.1 103.2 103.3]];  %rebuttal runs; DNA, DNA/Parb, DNA/Parb

expstotestrun=[102.1];  %rebuttal-2 runs; DNA, DNA/Parb, DNA/Parb

dryrun=0; %set to 1 and run only X050 to just find paths and files but spend no time on saving
%% start dipimage and run the experiment list
close all
X001_Dipstart;
disp('test running on:');
for ii=1:length(expstotestrun)
    expno=expstotestrun(ii);
    pathinfo = X000_setpath4snapshots(expno,usr);
    disp(pathinfo.dirmeasurement);
    %overwrite for test purposes:
    pathinfo.maxcellsperframe=2E10; %min 2 otherwise crash 
    %main programs in autrun test mode:
    if 1, X020_cellidentification4snapshots(pathinfo); end
    if 1, X030_cellcrop4snapshots(pathinfo); end
    if 1, X040_gridgeneration4snapshots(pathinfo); end
    if 1, X050_dnadensityuse(pathinfo,dryrun);    end
end
disp 'testrun ok for these datasets';
