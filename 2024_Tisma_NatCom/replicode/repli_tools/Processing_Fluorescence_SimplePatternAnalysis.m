function [fluo,modelpic]=Processing_Fluorescence_SimplePatternAnalysis(FL);
%JWJK_B:-------------------------------------------------------------------
%Title: Processing_Fluorescence_PatternAnalysis
%
%Approach: run it in autorun mode on provided test images to see.
%
%Input: FL_ori: image
%
%Output: structure 'fluo' with fields:
%     area_bac: 238
%     area_spot: 24
%     backbone_I: [1x31 double]
%     backbone_av: 4131.2
%     backbone_sigma: 1007.1
%     backbone_x: [1x31 double]
%     backbone_y: [1x31 double]
%     content_cytoplasm1: 3.5608e+005
%     content_signal: 4.1605e+005
%     content_spots1: 59972
%     content_total: 1426534
%     curve_medianofmax: [1x31 double]
%     curve_medianofmax_yposses: [1x31 double]
%     curve_medianofsum: [1x31 double]
%     level_dark: 931.13
%     level_edgetreshold: 1383.6
%     level_fluotreshold: 1157.4
%     level_medianofmax: 3726.6
%     level_medianofsum: 13312
%     level_peak: 7552
%     level_sumcyto: 13312
%     noise_dark: 113.11
%     peak_xpos: 6
%     peak_ypos: 17
%     ratio_FS: 0.14414
%     ratio_SN: 15.626
%     wherebac: [238x1 double]
%     wheredark: [631x1 double]
%     wherefluo: [361x1 double]
%     wherespot: [24x1 double]
%
% %References: M. Charl Moolman*, Jacob W.J. Kerssemakers*, and Nynke H. Dekker
% Quantitative analysis of intracellular fluorescent foci in live bacteria
% Biophysical Journal, online publication September 4 (2015)
%written by %JacobvKerssemakers, 2012
%
%:JWJK_B-------------------------------------------------------------------

    fluo.content_total=sum(FL(:));    
    fluo=Get_DarkProps(FL,fluo);   
    [fluo,modelpic]=Build_ModelPic(fluo,FL);
    fluo=orderfields(fluo);

function fluo=Get_DarkProps(FL,fluo);
    [r,c]=size(FL);
    %1 First, we make estimates on intensity and noise of the background (far outside
    %the bacterium).'local background' is defined as the average of the outer two image lines'
    if c>1
        fluo.level_dark=(mean(mean(FL(1:2,:)))+mean(mean(FL(r-1:r,:))))/2;
        %2) estimate dark (camera) noise via the spread in the difference between
        %neighbouring pixels
        diftop=FL(1:5,2:end)-FL(1:5,1:end-1); 
        difbot=FL(r-4:r,2:end)-FL(r-4:r,1:end-1);
        dif=[diftop(:);  difbot(:)];
        fluo.noise_dark=std(dif)/2^0.5;
        %define as 'fluorescence' those pixels sufficiently above the darklevel.
        %Note this may not be representative for the outline of the bacterium,
        %since there is some blurring and we want to measure all fluorescence
        fluotreshold=fluo.level_dark+2*fluo.noise_dark;
        %fluotreshold=fluo.level_dark-2*fluo.noise_dark;  %use this for montage making of DnaQ
        fluo.level_fluotreshold=fluotreshold;
        fluo.level_edgetreshold=fluo.level_dark+4*fluo.noise_dark;
    else
        fluo.noise_dark=0;
        fluo.level_dark=0;
        fluo.level_fluotreshold=0;
        fluo.level_edgetreshold=0;
    end

function fluo=Get_BackboneProps(FL, fluo);
    %first, we define a backbone of this bacterium:
    [r,c]=size(FL);
    if c>1
    bb0_x=[1:c]; axy=[1:c]; [X,Y]=meshgrid(bb0_x,axy);
    bb0_y=fluo.curve_medianofmax_yposses;
    bb0_I=fluo.curve_medianofmax;

    %select the higher part; build a picture with one backbone contour that
    %with an intensity of the median (y-summed) fluorescence counts; these
    %counts represent the cytoplasmic counts
    sel=find((bb0_I>fluo.level_fluotreshold-fluo.level_dark)...
              &(bb0_y>r/2-4)&(bb0_y<r/2+4)); 
          %includes xort-of centered
    lox=min(sel); hix=max(sel);
    fluo.backbone_x=bb0_x(lox:hix); 
    fluo.backbone_y=bb0_y(lox:hix); 
    fluo.backbone_I=bb0_I(lox:hix);    
    levels=0*fluo.backbone_x;
     for i=1:length(fluo.backbone_x)      
        %Get local area props!
        midx=fluo.backbone_x(i); midy=fluo.backbone_y(i);
        lox=max([midx-1,1]); hix=min([midx+1,c]);
        loy=max([midy-1,1]); hiy=min([midy+1,r]);
        squ=FL(loy:hiy,lox:hix)-fluo.level_dark;
        levels(i)=mean(squ(:)); 
     end
    fluo.backbone_av=mean(levels);
    fluo.backbone_sigma=std(levels);
    else
        fluo.backbone_x=0;
        fluo.backbone_y=0;
        fluo.backbone_I=0;
        fluo.backbone_av=0;
        fluo.backbone_sigma=0;
    end

function [fluo,modelpic]=Build_ModelPic(fluo,FL);
    %Build ModelPicmake a picture representing cytosol, and spot area, to show if the selected
    %area cover realistic parts; store these indices for later use
    modelpic=0*FL;
    [r,c]=size(FL);
        if c>1         
        fluo.wherefluo=find(FL>fluo.level_fluotreshold);
        fluo.wheredark=find(FL<=fluo.level_fluotreshold);
        fluo.wherebac=find(FL>fluo.level_edgetreshold);
        fluo.area_bac=length(fluo.wherebac);
        modelpic(fluo.wherefluo)=1;
        modelpic(fluo.wherebac)=2;
        else
        fluo.wherebac=[];
        fluo.wheredark=[];
        fluo.wherefluo=[];     
        fluo.area_bac=0;
        fluo.area_spot=0;
        end
    



    
   
