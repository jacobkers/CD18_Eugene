function [fluo,modelpic]=Processing_Fluorescence_PatternAnalysis(FL);
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


if nargin<1 %TEST MODUS;  
        close all
        pth=strcat(pwd, '\TestImages\');
        FL=1000+double(imread(strcat(pth,'TEST_SingleFocus.tif')));
 end       

    fluo.content_total=sum(FL(:));    
    fluo=Get_DarkProps(FL,fluo);
    fluo=Build_Fluorescence_vs_X_Curves(fluo,FL);   
    [fluo,modelpic]=Build_ModelPic(fluo,FL);
    fluo=Get_Peak_Props(fluo,FL);     
    fluo=Get_SpotsExcessContent(FL,fluo);
    fluo=Get_BackboneProps(FL, fluo);
    fluo=orderfields(fluo);


if  nargin<1
    fluo
    subplot(2,2,1); pcolor(FL); shading flat; 
    title('original image');
    subplot(2,2,2); pcolor(modelpic);
    title('regions: background-fluorescence-edge-spots');
    subplot(2,2,3); 
        plot(fluo.curve_medianofsum, 'k-o');hold on; 
        dummy=0*fluo.curve_medianofsum;
        plot(dummy+fluo. level_medianofsum,'r-');
        axis([0 length(sum(FL)) 0 1.2*max(fluo.curve_medianofsum)]);
        title('sums of X-sections');
        xlabel('x-position (pixel units)');
        ylabel('Summed Intensity')
    subplot(2,2,4); 
        plot(fluo.curve_medianofmax, 'k-o');hold on; 
        dummy=0*fluo.curve_medianofmax;
        plot(dummy+fluo. level_medianofmax,'r-');
        axis([0 length(sum(FL)) 0 1.2*max(fluo.curve_medianofmax)]);
        title('maxima of X-sections');
        xlabel('x-position (pixel units)');
        ylabel('Max.Intensity')
    [~]=ginput(1);
    close(gcf);
end
end %if test


 function fluo=Get_SpotsExcessContent(FL,fluo);
        [r,c]=size(FL);
        if c>1
            [X,Y]=meshgrid(1:c,1:r); 
            bacminxpos=min(X(fluo.wherebac));
            bacmaxpos=max(X(fluo.wherebac));
            fluo.level_sumcyto=median(fluo.curve_medianofsum(bacminxpos:bacmaxpos));
            sel=find(fluo.curve_medianofsum>fluo.level_sumcyto);
            fluo.content_spots1=sum(fluo.curve_medianofsum(sel)-fluo.level_sumcyto);
            fluo.content_cytoplasm1=fluo.content_signal-fluo.content_spots1;
            %Last, detrmine some handy ratios
            fluo.ratio_FS=fluo.content_spots1/fluo.content_signal;
            fluo.ratio_SN=mean(fluo.curve_medianofmax)/(2*fluo.noise_dark);
        else
            fluo.content_cytoplasm1=0;
            fluo.content_signal=0;
            fluo.content_spots1=0;
            fluo.content_total=0;
            fluo.ratio_FS=0;
            fluo.ratio_SN=0; 
        end
     end

function fluo=Get_Peak_Props(fluo,FL);
    [r,c]=size(FL);
    if c>1
        %find global max as crude spot location
        fluo.level_peak=max(FL(:));
        [~,fluo.peak_xpos]=max(max(FL));
        [~,fluo.peak_ypos]=max(max(FL'));
    else 
         fluo.level_peak=0;
        fluo.peak_xpos=0;
        fluo.peak_ypos=0;   
    end
end


function fluo=Build_Fluorescence_vs_X_Curves(fluo,FL);
    %make two 1D-curves along the x direction (length of bacterium): 
    %-one that contains all the counts; 
    %one that contains the maxima along the y-direction.
    %these curves are used to find 'median' levels, representative for the
    %cytoplasmic signal level
    [r,c]=size(FL);
    if c>1
    backsel=find(FL<fluo.level_fluotreshold);
    FLbc=FL-fluo.level_fluotreshold;
    FLbc(backsel)=0;    
    [X,~]=meshgrid(1:c,1:r);
    
    %a) First, the summed intensities (minus background)
    fluo.curve_medianofsum=sum(FLbc);  %y-summed,including out-of bacteria-edges
    fluo.level_medianofsum=median(fluo.curve_medianofsum);
    fluo.content_signal=sum(fluo.curve_medianofsum);       %total represents total counts
        
    %b) As a second approach to counting spot content, we use the 'maxima' median
    %level and the full picture. This 2D-approach may give us a sharper distinction between
    %'spot' and 'no-spot' area
    [fluomaxcurve,maxidx]=max(FLbc);
    fluo.curve_medianofmax=fluomaxcurve;
    fluo.level_medianofmax=median(fluomaxcurve);  %level represents backbone level bacterium
    fluo.curve_medianofmax_yposses=maxidx;  %line showing position of maxima, spanning whole picture
    else
    fluo.curve_medianofmax=1;
    fluo.curve_medianofsum=1;
    fluo.level_medianofmax=0;
    fluo.level_medianofmax_yposcurve=0;
    fluo.level_medianofsum=0;
    end
end 


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
end

function [fluo,modelpic]=Build_ModelPic(fluo,FL);
    %Build ModelPicmake a picture representing cytosol, and spot area, to show if the selected
    %area cover realistic parts; store these indices for later use
    modelpic=0*FL;
    [r,c]=size(FL);
        if c>1         
        fluo.wherefluo=find(FL>fluo.level_fluotreshold);
        fluo.wheredark=find(FL<=fluo.level_fluotreshold);
        fluo.wherespot=find(FL>fluo.level_medianofmax+fluo.level_fluotreshold);
        fluo.wherebac=find(FL>fluo.level_edgetreshold);
        fluo.area_bac=length(fluo.wherebac);
        fluo.area_spot=length(fluo.wherespot);  %amount of 'spot pixels'. 
        modelpic(fluo.wherefluo)=1;
        modelpic(fluo.wherebac)=2;
        modelpic(fluo.wherespot)=3;
        else
        fluo.wherebac=[];
        fluo.wheredark=[];
        fluo.wherefluo=[];
        fluo.wherespot=[];     
        fluo.area_bac=0;
        fluo.area_spot=0;
        end
end
    



    
   
