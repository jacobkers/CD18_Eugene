function [info_Cnd_allROIs,info_DNA_allROIs]=spots00_harvest_all_rois(expi,init,AllExp,outpath);
%JWJK_B:----[add ABCorC*----------------------------------------------------
%Title: Gather all info of different tethers
%Summary: Collect all spot position ans content data from all chosen rois, 
%to make histograms etc.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_B-----[add ABCorC*---------------------------------------------------
%
info_Cnd_allROIs=struct(...
        'pos_roino',[],...
        'pos_frameno',[],...
        'pos_X_pix',[],...
        'pos_X_subpix',[],...
        'content_peakvals',[],...
        'content_perspot_est',[],...
        'content_perspot_meas',[],...
        'free_spots_numbers',[],...
        'free_tether_lengths',[]...
        );
    
 info_DNA_allROIs=struct(...
        'pos_roino',[],...
        'pos_frameno',[],...
        'pos_X_pix',[],...
        'pos_X_subpix',[],...
        'content_peakvals',[],...
        'content_perspot_est',[],...
        'content_perspot_meas',[]...
        );
    
LE=length(AllExp);  %for all experiments
SaveName=char(strcat(outpath, 'AllROI_allresults.mat'));
 close all;
for roi=1:LE  
    Exp=strcat(init.roidirname,num2str(AllExp(roi)));
    expinfo=A002_Condensin_with_plectonemes_expinfo(expi,AllExp(roi));
    corr=expinfo.channelshift;    
    LoadName=char(strcat(outpath, Exp));             
    disp(strcat('Harvesting data: Exps to work through:',num2str(LE-roi)));
    load(strcat(LoadName, '_allresults.mat')); 
    
    Ld=length(info_DNA.pos_frameno);
    info_DNA.pos_roino=1:Ld;
    
    info_DNA_allROIs.pos_roino=[info_DNA_allROIs.pos_roino info_DNA.pos_roino];
    info_DNA_allROIs.pos_frameno=[info_DNA_allROIs.pos_frameno info_DNA.pos_frameno];
    info_DNA_allROIs.pos_X_pix=[info_DNA_allROIs.pos_X_pix info_DNA.pos_X_pix];
    info_DNA_allROIs.pos_X_subpix=[info_DNA_allROIs.pos_X_subpix info_DNA.pos_X_subpix];
    info_DNA_allROIs.content_peakvals=[info_DNA_allROIs.content_peakvals info_DNA.content_peakvals];
    info_DNA_allROIs.content_perspot_est=[info_DNA_allROIs.content_perspot_est info_DNA.content_perspot_est];
    info_DNA_allROIs.content_perspot_meas=[info_DNA_allROIs.content_perspot_meas info_DNA.content_perspot_meas];
    
    Lc=length(info_Cnd.pos_frameno);
    info_Cnd.pos_roino=1:Lc;
    info_Cnd_allROIs.pos_roino=[info_Cnd_allROIs.pos_roino info_Cnd.pos_roino];
    info_Cnd_allROIs.pos_frameno=[info_Cnd_allROIs.pos_frameno info_Cnd.pos_frameno];
    info_Cnd_allROIs.pos_X_pix=[info_Cnd_allROIs.pos_X_pix info_Cnd.pos_X_pix+corr];
    info_Cnd_allROIs.pos_X_subpix=[info_Cnd_allROIs.pos_X_subpix info_Cnd.pos_X_subpix+corr];
    info_Cnd_allROIs.content_peakvals=[info_Cnd_allROIs.content_peakvals info_Cnd.content_peakvals];
    info_Cnd_allROIs.content_perspot_est=[info_Cnd_allROIs.content_perspot_est info_Cnd.content_perspot_est];
    info_Cnd_allROIs.content_perspot_meas=[info_Cnd_allROIs.content_perspot_meas info_Cnd.content_perspot_meas];
    info_Cnd_allROIs.free_spots_numbers=[info_Cnd_allROIs.free_spots_numbers  info_Cnd.general_free_number];
    info_Cnd_allROIs.free_tether_lengths=[info_Cnd_allROIs.free_tether_lengths info_Cnd.general_total_freetetherlength];
end
info_Cnd_allROIs.free_spots_densities=info_Cnd_allROIs.free_spots_numbers./...
                                      info_Cnd_allROIs.free_tether_lengths;
info_Cnd_allROIs.free_spots_densities_av=nanmean(info_Cnd_allROIs.free_spots_densities);
info_Cnd_allROIs.free_spots_densities_st=nanstd(info_Cnd_allROIs.free_spots_densities);

dum=1;
save(SaveName,'info_DNA_allROIs','info_Cnd_allROIs'); 
