function kymo_1=kym_convert_dna_kymo(kymo_DNA,levels_DNA);
%JWJK_C:----[add ABCorC*----------------------------------------------------
%Title: %convert fluorescene kymograph to genomic percentage. 
%Summary: Simple normalization routine, counting all fluorescence counts
%above a certain noise treshold. All counts define 100% genome content.
%Input: kymograph.
%Output: normalized kymograph.
%References: CD lab, project Eugene Kim, written by Jacob Kers, 2019
%:JWJK_C-----[add ABCorC*---------------------------------------------------


% levels_DNA
%     level_dark: 1.2870
%     noise_dark: 0.3811
%     tetherlevel: 6.9567
%     level_darktreshold: 2.0493
%     level_tethertreshold: 6.1945
%     level_looptreshold: 7.7190

kymo_1=kymo_DNA-levels_DNA.level_darktreshold;
kymo_1(kymo_1<0)=0;
[rr,cc]=size(kymo_1);
for ii=1:rr
    kymo_1(ii,:)=kymo_1(ii,:)/sum(kymo_1(ii,:))*100;
end
             
if 0
    subplot(1,2,1); 
        pcolor(kymo_DNA); shading flat; colormap bone; hold on;
        title('DNA');
    subplot(1,2,2); 
        pcolor(kymo_1); shading flat; colormap bone; hold on;
        title('perc');
    dum=1;
end