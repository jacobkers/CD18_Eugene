function plot_content_per_frame
%this stub loads processed data and plots content per frame
clc; close all;


%% Setting up the user paths
override_userpaths_to_sharedpaths=1;  
%setting this to one overwrites local paths and uses shared paths (be sure
%the raw data is there)
switch 1
    case 1, usr='Eugene';  
        expi=9; %'2020_01_13 MukBEF_msd_wATP\'
        init=A001_EK_condensin_with_plectonemes_init(expi);
    case 2, usr='Jacob';
        expi=0; %'2019_09_02_NumberOfCondensinPerDNA\'  TEST 
        expi=3; %'2019_09_02_NumberOfCondensinPerDNA\'  all rois 
        init=A001_JK_condensin_with_plectonemes_init(expi);
end

if override_userpaths_to_sharedpaths
   mpth='M:\tnw\bn\cd\Shared\Jacob Kerssemakers\';  
   init.datapathin=[mpth 'TESTdata_out\Eugene\matlabresults\'];          %loading
   init.datapathout=[mpth 'TESTdata_out\Eugene\matlabresults\']; %saving
   addpath(genpath(swap_path('Dropbox\CD_recent\BN_CD18_Eugene\Matlab\Matlab_tools_CD18EK\'))); %tools:
end

savepath=[init.datapathout init.expname, 'contentplots\'];
% if ~isdir(savepath), mkdir(savepath); end

%% now, load data
pth=[init.datapathin init.expname, 'collected\'];
load([pth,'EKMcp_A040_AllROI_plectoneme_free.mat']);
load([pth, 'EKMcp_A040_AllROI_plectoneme_Cnd_associated.mat']);

N_rois=length(init.AllExp);

for ii=1:N_rois
roino=init.AllExp(ii);

%get content per frame
[content_plec_free,N_free]=get_content_per_frame(info_DNA_free,roino);
[content_plec_cnd,N_cnd]=get_content_per_frame(info_DNA_near_Cnd,roino);

 savename=['content_roi', num2str(roino, '% 02.0f'),'.jpg']; 

subplot(2,1,1);
    plot(content_plec_free); hold on;
    plot(content_plec_cnd, 'r');
    title(Replace_underscores(savename(1:end-4)));
    legend('free', 'with-condensin');
    xlabel('frame');
    ylabel('perc');
subplot(2,1,2);
    plot(N_free); hold on;
    plot(N_cnd, 'r');
    legend('free', 'with-condensin');
    xlabel('frame');
    ylabel('# plectonemes');

 saveas(gcf,[savepath savename]);
 pause(0.1);
 close(gcf);
end

function [content, no_of_peaks]=get_content_per_frame(info,roino);
%get selected content per frame
N_roi=max(info.pos_roino);

%pick one roi

roi_sel=find(info.pos_roino==roino);
pks=info.content_peakvals(roi_sel);  %peak content
frs=info.pos_frameno(roi_sel);         %frames

N_frames=max(frs);
content=NaN*zeros(N_frames,1);
no_of_peaks=zeros(N_frames,1);
for fi=1:N_frames
    fr_sel=find(frs==fi);
    if ~isempty (fr_sel)
        content(fi)=sum(pks(fr_sel)); %sum peaks per frame
        no_of_peaks(fi)=length(fr_sel); %# peaks per frame
    end
end

