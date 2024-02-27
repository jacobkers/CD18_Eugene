function histSpotsNumber_jk_adapt_bootstrap
%% Written by Géraldine Laloux, UCLouvain, Feb 2020
%% Microscope resolution (um to pixels conversion) - modify this value when required
%px2um = 0.07;  Dumbledore 1.5X pixel size
%px2um = 0.11;  Dumbledore pixel size
%px2um = 0.04;  Snape 1.5X pixel size
%px2um = 0.065; Snape pixel size
clc; close all;
cur_path=[pwd,'\'];
switch 4
    case 1 %local test 1
        source_path=pwd;
        nme='4595_30min_002c1_Spotdetection';
    case 2 %local test 2
        source_path=pwd;
        nme='4595_180min_001c2_Spotdetection';
    case 3 %16-10-23 figure 2c rebuttal
        source_path='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\20220706_4595_kale_1\'
        nme='4595_IPTG_0min_c1_Stack_Corrected_CellOutlines_SpotDetection_01_final';
    case 4 %16-10-23 figure 2c rebuttal
        source_path='M:\tnw\bn\cd\Shared\Tisma\Bacillus subtilis chromosome imaging\Widefield_and_deconvolution\20221004_IPTG_halt_timing\'
        nme='4595_180min_001c2_Spotdetection';
end

load([source_path nme '.mat']);

px2um = 0.065; 
%% This section gets extradata from the meshData provided by Oufti (to obtain calculated length, width, ... among other cell features)
% and creates a new cellList named cellListExtra
for frame=1:length(cellList.meshData)
    for cell=1:length(cellList.meshData{frame})
        if isfield (cellList.meshData{frame}{cell},'length')
            cellListExtra.meshData{frame}{cell} = cellList.meshData{frame}{cell};
        end
        if ~isfield (cellList.meshData{frame}{cell},'length')
        cellListExtra.meshData{frame}{cell} = getextradata(cellList.meshData{frame}{cell});
        end
    end
end
%% Script
n = 0;
cellL = cellListExtra.meshData;
spotNumber = [];
cellLength = [];
for frame=1:length(cellL)
    for cell=1:length(cellL{frame})
        if isfield(cellL{frame}{cell},'spots')
           
            n = n + 1;
      
            spots = length(cellL{frame}{cell}.spots.l);
            spotNumber = [spotNumber spots];
            
        end
    end
end


%% Histogram settings and bootstrapping
bin_ax = 0:1:8; % this gives the range and step size of the histogram (ex: 0:1:15 = from 0 to 15 with steps of 1). Modify it if needed.
hist_counts = hist(spotNumber,bin_ax);
[bootstat,bootsam_idx] = bootstrp(1000,@mean,spotNumber);
boot_samples=spotNumber(bootsam_idx);
h_Nboot = hist(boot_samples,bin_ax);
hist_perc = 100*hist_counts/(sum(hist_counts));
hPERC_err=2* 100*std(h_Nboot')/sum(hist_counts);  %95%
%% Figure
figure1 = figure;
bar(bin_ax,hist_perc,'k'); hold on; % this gives the range and step size of the histogram (ex: 0:1:15 = from 0 to 15 with steps of 1). Modify it if needed.
errorbar(bin_ax,hist_perc,hPERC_err, 'rsq', 'LineWidth', 2);
text(6, 50,['N=' num2str(n)]);
xlabel('Number of foci per cell','FontSize',20);
ylabel('Fraction of cells (%)','FontSize',20);
box off
ylim ([0 100]); % Limits of Y axis. Change if needed.
%% Display number of cells 
disp(['cell number:' num2str(n)]);
disp(['mean spots: ' num2str(mean(spotNumber)) ' ± ' num2str(std(spotNumber))])
%% Saving as pdf , fig , tables
outdirs=[{cur_path}, {source_path}];
for dd=1:2
    outdir=[outdirs{dd} , '\bootstrapped\'];
    mkdir(outdir);
    saveas(gcf,[outdir nme '_nspots.jpg']);
    saveas(gcf,[outdir nme '_nspots.pdf']);
    saveas(gcf,[outdir nme '_nspots.fig']);
    table_hdrs=[{'bins'}, {'counts'}, {'percentage'},{'error'}];
    table_data=[bin_ax',hist_counts', hist_perc',hPERC_err'];
    cd(outdir);
    %xlswrite(['Wielontwerpen\', ontwerp, '_maten.xlsx'], ColNames,'snijplan per strook','A1');
    xlswrite([nme '_nspots.xlsx'],table_hdrs,'bootstrapped','A1');
    xlswrite([nme '_nspots.xlsx'],table_data,'bootstrapped','A2');
    cd(cur_path);
end
%tqbualr results