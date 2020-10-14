function [info_Cnd_allROIs,info_DNA_allROIs,info_Cnd_per_ROI,info_DNA_per_ROI]=spots00_harvest_all_ROIs(expi,init,AllExp,usr);
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
        'content_clustercont',[],...
        'content_perspot_est',[],...
        'content_perspot_meas',[],...
        'free_spots_numbers',[],...
        'free_tether_lengths',[]);
 info_Cnd_allROIs.label.farfrom_dna_edges=[];
    
 info_DNA_allROIs=struct(...
        'pos_roino',[],...
        'pos_frameno',[],...
        'pos_X_pix',[],...
        'pos_X_subpix',[],...
        'content_clustercont',[],...
        'content_perspot_est',[],...
        'content_perspot_meas',[]);
 info_DNA_allROIs.label.farfrom_dna_edges=[];
    
    
LE=length(AllExp);  %for all experiments
inpath=strcat(init.datapathout, 'matlabresults\',init.expname,'\');

close all;
for roi=1:LE  
    Exp=strcat(init.roidirname,num2str(AllExp(roi)));
    roino=AllExp(roi);
    switch usr
    case 'Jacob',  expinfo=A002_JK_Condensin_with_plectonemes_expinfo(expi,AllExp(roi));
    case 'Eugene', expinfo=A002_EK_Condensin_with_plectonemes_expinfo(expi,AllExp(roi));
    end    
    
    corr=expinfo.channelshift;    
    LoadName=char(strcat(inpath, 'EKMcp_A020_',Exp));             
    disp(strcat('Harvesting data: Exps to work through:',num2str(LE-roi)));
    load(strcat(LoadName, '_allresults.mat'));     
    
    
    [kymo_duration,kymo_width]=size(kymo_DNA);
     
    info_DNA_per_ROI.kymo_duration(roi)=kymo_duration;
    info_DNA_per_ROI.kymo_width(roi)=kymo_width;
    info_DNA_per_ROI.channelshift(roi)=expinfo.channelshift;
    info_DNA_per_ROI.SaveName{roi}=char(strcat('Roi_',Exp)); 
    
    info_Cnd_per_ROI.kymo_duration(roi)=kymo_duration;
    info_Cnd_per_ROI.kymo_width(roi)=kymo_width;  
    info_Cnd_per_ROI.channelshift(roi)=expinfo.channelshift;
    info_Cnd_per_ROI.SaveName{roi}=char(strcat('Roi_',Exp)); 
    
    Ld=length(info_DNA.pos_frameno);
    pos_roino_DNA=0*(1:Ld)+roino; 
    info_DNA_allROIs.pos_roino=[info_DNA_allROIs.pos_roino pos_roino_DNA];
    info_DNA_allROIs.pos_frameno=[info_DNA_allROIs.pos_frameno info_DNA.pos_frameno];
    info_DNA_allROIs.pos_X_pix=[info_DNA_allROIs.pos_X_pix info_DNA.pos_X_pix];
    info_DNA_allROIs.pos_X_subpix=[info_DNA_allROIs.pos_X_subpix info_DNA.pos_X_subpix];
    info_DNA_allROIs.content_clustercont=[info_DNA_allROIs.content_clustercont info_DNA.content_clustercont];
    info_DNA_allROIs.content_perspot_est=[info_DNA_allROIs.content_perspot_est info_DNA.content_perspot_est];
    info_DNA_allROIs.content_perspot_meas=[info_DNA_allROIs.content_perspot_meas info_DNA.content_perspot_meas];
    info_DNA_allROIs.label.farfrom_dna_edges=[info_DNA_allROIs.label.farfrom_dna_edges info_DNA.classify.awayfromDNAedges];
    
    Lc=length(info_Cnd.pos_frameno);
    pos_roino_Cnd=0*(1:Lc)+roino;Lc=length(info_Cnd.pos_frameno); 
    info_Cnd_allROIs.pos_roino=[info_Cnd_allROIs.pos_roino pos_roino_Cnd];
    info_Cnd_allROIs.pos_frameno=[info_Cnd_allROIs.pos_frameno info_Cnd.pos_frameno];
    info_Cnd_allROIs.pos_X_pix=[info_Cnd_allROIs.pos_X_pix info_Cnd.pos_X_pix+corr];
    info_Cnd_allROIs.pos_X_subpix=[info_Cnd_allROIs.pos_X_subpix info_Cnd.pos_X_subpix+corr];
    info_Cnd_allROIs.content_clustercont=[info_Cnd_allROIs.content_clustercont info_Cnd.content_clustercont];
    info_Cnd_allROIs.content_perspot_est=[info_Cnd_allROIs.content_perspot_est info_Cnd.content_perspot_est];
    info_Cnd_allROIs.content_perspot_meas=[info_Cnd_allROIs.content_perspot_meas info_Cnd.content_perspot_meas];
    info_Cnd_allROIs.free_spots_numbers=[info_Cnd_allROIs.free_spots_numbers  info_Cnd.general_free_number];
    info_Cnd_allROIs.free_tether_lengths=[info_Cnd_allROIs.free_tether_lengths info_Cnd.general_total_freetetherlength];
    info_Cnd_allROIs.label.farfrom_dna_edges=[info_Cnd_allROIs.label.farfrom_dna_edges info_Cnd.classify.awayfromDNAedges];
end
info_Cnd_allROIs.free_spots_densities=info_Cnd_allROIs.free_spots_numbers./...
                                      info_Cnd_allROIs.free_tether_lengths;
info_Cnd_allROIs.free_spots_densities_av=nanmean(info_Cnd_allROIs.free_spots_densities);
info_Cnd_allROIs.free_spots_densities_st=nanstd(info_Cnd_allROIs.free_spots_densities);

dum=1;

